"""
HOTPANTS High-Level Python API

Provides user-friendly dataclasses and functions for kernel fitting and
image differencing. This module orchestrates the low-level cffi calls
and handles numpy array I/O.

Example usage:

    import numpy as np
    from hotpants import fit_kernel, spatial_convolve

    # Load images
    template = np.load('template.npy').astype(np.float32)
    science = np.load('science.npy').astype(np.float32)

    # Configure and fit kernel
    config = KernelConfig()
    layout = RegionLayout()
    thresholds = NoiseThresholds()

    solution = fit_kernel(template, science, noise=None,
                          config=config, layout=layout, thresholds=thresholds)

    # Apply kernel to create difference image
    difference = spatial_convolve(science, solution, config)

Reference: Alard & Lupton (1998), ApJ 503:325
"""

import numpy as np
from loguru import logger
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from . import _core, _hotpants_ffi
from ._core import allocate_array, global_state, validate_image_array

# =====================================================================
# Configuration Models (Pydantic)
# =====================================================================


class KernelConfig(BaseModel):
    """
    Kernel fitting and basis configuration.

    Attributes:
        kernel_half_width: Half-width of kernel region (pixels). Must be positive.
            Controls the spatial scale of the kernel. Typical: 10-20 px.
            (Maps to C global: hwKernel)

        kernel_order: Spatial polynomial order for kernel coefficients.
            Controls how the kernel varies with position. Typical: 2-4.
            Higher values fit more spatial variation but risk overfitting.
            (Maps to C global: kerOrder)

        bg_order: Spatial polynomial order for background (sky) level.
            Typical: 1-2. Set to 0 for constant background.
            (Maps to C global: bgOrder)

        fit_threshold: Sigma threshold for stamps used in kernel fit.
            Only pixels > mean + fit_threshold * std are used. Typical: 20.
            (Maps to C global: kerFitThresh)

        scale_fit_threshold: Scale factor applied to fit_threshold in retries.
            When initial fit fails, reduced to fit_threshold * scale_fit_threshold.
            Typical: 0.5.
            (Maps to C global: scaleFitThresh)

        n_ks_stamps: Number of kernel test substamps per stamp.
            Bright stars in each stamp provide fit constraints. Typical: 3.
            (Maps to C global: nKSStamps)

        hw_ks_stamp: Half-width of kernel test substamps (pixels).
            Must be <= kernel_half_width. Typical: 10.
            (Maps to C global: hwKSStamp)
    """

    kernel_half_width: int = Field(default=15, gt=0)
    kernel_order: int = Field(default=2, ge=0)
    bg_order: int = Field(default=1, ge=0)
    fit_threshold: float = Field(default=20.0, gt=0)
    scale_fit_threshold: float = Field(default=0.5, gt=0, le=1)
    n_ks_stamps: int = Field(default=3, gt=0)
    hw_ks_stamp: int = Field(default=10, gt=0)

    @field_validator("hw_ks_stamp")
    @classmethod
    def validate_ks_stamp_size(cls, v: int, info) -> int:
        """hw_ks_stamp is independent of kernel_half_width; just ensure positive."""
        return v


class RegionLayout(BaseModel):
    """
    Image tiling and stamp grid configuration.

    Attributes:
        n_regions_x: Number of regions to divide image into (X axis).
            Each region is fit separately, allowing for spatially-varying kernel.
            Typical: 1-3. Larger values for wide-field images.
            (Maps to C global: nRegX)

        n_regions_y: Number of regions (Y axis). Usually matches n_regions_x.
            (Maps to C global: nRegY)

        stamps_per_region_x: Number of stamp grid points per region (X axis).
            Each stamp is a small cutout where the kernel is fit. Typical: 10.
            (Maps to C global: nStampX)

        stamps_per_region_y: Number of stamps per region (Y axis).
            (Maps to C global: nStampY)

        use_full_substamps: If True, automatically divide region into kernel-sized
            substamps (ignores stamps_per_region_*). Typical: False.
            (Maps to C global: useFullSS)
    """

    n_regions_x: int = Field(default=1, ge=1)
    n_regions_y: int = Field(default=1, ge=1)
    stamps_per_region_x: int = Field(default=10, ge=1)
    stamps_per_region_y: int = Field(default=10, ge=1)
    use_full_substamps: bool = False


class NoiseThresholds(BaseModel):
    """
    Data quality thresholds and noise model.

    Attributes:
        template_upper_threshold: Upper valid data threshold for template.
            Pixels > this are considered saturated and masked. Typical: 25000.
            (Maps to C global: tUThresh)

        template_lower_threshold: Lower valid data threshold for template.
            Pixels < this are masked. Typical: 0.
            (Maps to C global: tLThresh)

        science_upper_threshold: Upper valid threshold for science image.
            (Maps to C global: iUThresh)

        science_lower_threshold: Lower valid threshold for science image.
            (Maps to C global: iLThresh)

        template_gain: Gain of template (e-/ADU). Used for noise estimation.
            Typical: 1.0-3.0. Must be positive.
            (Maps to C global: tGain)

        science_gain: Gain of science image (e-/ADU).
            (Maps to C global: iGain)

        template_readnoise: Template read noise (electrons). Typical: 0-100 e-.
            (Maps to C global: tRdnoise)

        science_readnoise: Science image read noise (e-).
            (Maps to C global: iRdnoise)

        template_pedestal: Baseline offset (ADU) in template. Typical: 0.
            (Maps to C global: tPedestal)

        science_pedestal: Baseline offset in science image.
            (Maps to C global: iPedestal)
    """

    template_upper_threshold: float = 25000.0
    template_lower_threshold: float = 0.0
    science_upper_threshold: float = 25000.0
    science_lower_threshold: float = 0.0
    template_gain: float = Field(default=1.0, gt=0)
    science_gain: float = Field(default=1.0, gt=0)
    template_readnoise: float = Field(default=0.0, ge=0)
    science_readnoise: float = Field(default=0.0, ge=0)
    template_pedestal: float = 0.0
    science_pedestal: float = 0.0

    @model_validator(mode="after")
    def validate_thresholds(self) -> "NoiseThresholds":
        """Ensure upper thresholds > lower thresholds."""
        if self.template_upper_threshold <= self.template_lower_threshold:
            msg = "template_upper_threshold must be > template_lower_threshold"
            raise ValueError(msg)
        if self.science_upper_threshold <= self.science_lower_threshold:
            msg = "science_upper_threshold must be > science_lower_threshold"
            raise ValueError(msg)
        return self


class KernelSolution(BaseModel):
    """
    Result of kernel fitting.

    Attributes:
        chi2: Fit residual (chi-squared).
        kernel_norm: Kernel integral (sum).
        mean_sigma: Mean significance of stamps used in fit.
        scatter_sigma: Scatter in stamp significance.
        n_skipped_stamps: Number of stamps excluded due to poor fit.
        kernel_coefficients: Array of fitted kernel polynomial coefficients.
            Shape: (n_kernel_components,)
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    chi2: float
    kernel_norm: float
    mean_sigma: float
    scatter_sigma: float
    n_skipped_stamps: int
    kernel_coefficients: np.ndarray

    def __repr__(self) -> str:
        return (
            f"KernelSolution(chi2={self.chi2:.3f}, "
            f"norm={self.kernel_norm:.3f}, "
            f"mean_sigma={self.mean_sigma:.2f}, "
            f"n_skipped={self.n_skipped_stamps})"
        )


# =====================================================================
# High-Level API Functions
# =====================================================================


def fit_kernel(
    template: np.ndarray,
    science: np.ndarray,
    noise: np.ndarray | None = None,
    config: KernelConfig | None = None,
    layout: RegionLayout | None = None,
    thresholds: NoiseThresholds | None = None,
    verbose: int = 0,
    n_thread: int = 1,
) -> KernelSolution:
    """
    Fit a spatially-varying convolution kernel to match template and science PSFs.

    This function orchestrates the kernel fitting algorithm from Alard & Lupton (1998):
    1. Divide image into regions and lay out stamp grid
    2. Identify bright star centers in each stamp
    3. Accumulate normal equations for least-squares kernel fit
    4. Solve for spatially-varying kernel coefficients via Cholesky decomposition
    5. Return fitted kernel with quality metrics

    Args:
        template: Reference image (numpy array, float32, 2D).
            Should be well-registered and have good PSF stability.

        science: Science image to be differenced (numpy array, float32, 2D).
            Must have same shape as template.

        noise: Optional noise image (float32, 2D) or None.
            If provided, used to weight the fit. Typically variance or per-pixel RMS.
            If None, noise is estimated from data statistics.

        config: KernelConfig instance (default: KernelConfig()).
            Controls kernel order, fit thresholds, and basis configuration.

        layout: RegionLayout instance (default: RegionLayout()).
            Controls image tiling and stamp grid.

        thresholds: NoiseThresholds instance (default: NoiseThresholds()).
            Specifies saturation thresholds and noise model parameters.

        verbose: Verbosity level (0=silent, 1=progress, 2=debug).

        n_thread: Number of threads for parallelism (default: 1).
            Note: Requires OpenMP support in C library.

    Returns:
        KernelSolution: fitted kernel and quality metrics

    Raises:
        TypeError: if image arrays are not float32
        ValueError: if images don't match size, contain NaN/Inf, or
                    configuration parameters are invalid
        RuntimeError: if kernel fitting fails (e.g., singular matrix)

    Example:
        template = np.random.randn(512, 512).astype(np.float32)
        science = np.random.randn(512, 512).astype(np.float32)

        solution = fit_kernel(template, science)
        print(f"Fitted kernel chi2: {solution.chi2:.3f}")
        print(f"Kernel integral: {solution.kernel_norm:.3f}")
    """
    # Use defaults if not provided
    if config is None:
        config = KernelConfig()
    if layout is None:
        layout = RegionLayout()
    if thresholds is None:
        thresholds = NoiseThresholds()

    # Validate inputs
    validate_image_array(template, "template")
    validate_image_array(science, "science")

    if template.shape != science.shape:
        msg = f"Image shape mismatch: template {template.shape} != science {science.shape}"
        raise ValueError(msg)

    if noise is not None:
        validate_image_array(noise, "noise")
        if noise.shape != template.shape:
            msg = f"Noise shape mismatch: {noise.shape} != {template.shape}"
            raise ValueError(msg)

    ny, nx = template.shape

    logger.debug(f"Template shape: {ny},{nx}")

    # Convert to C-contiguous if needed
    if not template.flags["C_CONTIGUOUS"]:
        template = np.ascontiguousarray(template)
    if not science.flags["C_CONTIGUOUS"]:
        science = np.ascontiguousarray(science)
    if noise is not None and not noise.flags["C_CONTIGUOUS"]:
        noise = np.ascontiguousarray(noise)

    # Set global configuration
    config_dict = {
        "hwKernel": config.kernel_half_width,
        "kerOrder": config.kernel_order,
        "bgOrder": config.bg_order,
        "kerFitThresh": config.fit_threshold,
        "scaleFitThresh": config.scale_fit_threshold,
        "nKSStamps": config.n_ks_stamps,
        "hwKSStamp": config.hw_ks_stamp,
        "nRegX": layout.n_regions_x,
        "nRegY": layout.n_regions_y,
        "nStampX": layout.stamps_per_region_x,
        "nStampY": layout.stamps_per_region_y,
        "useFullSS": 1 if layout.use_full_substamps else 0,
        "tUThresh": thresholds.template_upper_threshold,
        "tLThresh": thresholds.template_lower_threshold,
        "iUThresh": thresholds.science_upper_threshold,
        "iLThresh": thresholds.science_lower_threshold,
        "tGain": thresholds.template_gain,
        "iGain": thresholds.science_gain,
        "tRdnoise": thresholds.template_readnoise,
        "iRdnoise": thresholds.science_readnoise,
        "tPedestal": thresholds.template_pedestal,
        "iPedestal": thresholds.science_pedestal,
        "verbose": verbose,
        "nThread": n_thread,
        "useTPS": 0,  # Will be set after global_state context is set up
        "tpsSmoothing": 0.0,
    }

    with global_state(config_dict):
        # Build stamp grid and compute statistics
        logger.debug("Building stamp grid...")
        stamps = _core.build_stamps(
            template,
            science,
            layout.n_regions_x,
            layout.n_regions_y,
            layout.stamps_per_region_x,
            layout.stamps_per_region_y,
        )
        # Note: build_stamps does NOT call cleanupBuildStampsContext().
        # The region context (mRData, region buffers) remains active because
        # fitKernel → check_again → getStampSig/fillStamp access mRData.
        # We must call cleanup after fitKernel completes (see finally block below).

        try:
            # Fit kernel
            logger.debug("Fitting kernel...")
            kernel_info = _core.get_kernel_info()
            # Calculate kernelSol size for polynomial or TPS mode
            # Note: TPS uses more space to store RBF weights and control point positions
            use_tps = False  # Currently TPS is not exposed to user; default to polynomial
            n_coeffs = _core.calculate_kernel_solution_size(
                layout.stamps_per_region_x,
                layout.stamps_per_region_y,
                use_tps=use_tps,
            )

            # If no noise image provided, compute Poisson variance to match main.c's
            # makeNoiseImage4 call: var = |img| / gain + (rdnoise / gain)^2.
            # A zero-filled array causes division-by-zero in getStampSig → all stamps
            # get rejected. Proper variance is required for stamp quality assessment.
            if noise is None:
                t_var = (
                    np.abs(template.astype(np.float64)) / thresholds.template_gain
                    + (thresholds.template_readnoise / thresholds.template_gain) ** 2
                )
                i_var = (
                    np.abs(science.astype(np.float64)) / thresholds.science_gain
                    + (thresholds.science_readnoise / thresholds.science_gain) ** 2
                )
                noise = (t_var + i_var).astype(np.float32)
            else:
                noise = np.ascontiguousarray(noise.astype(np.float32))

            (
                chi2,
                kernel_norm,
                mean_sigma,
                scatter_sigma,
                n_skipped,
                kernel_coeffs,
            ) = _core.fit_kernel_c(
                stamps,
                template,
                science,
                noise,
                n_coeffs,
            )

            logger.debug(f"Fitted kernel: chi2={chi2:.3f}, norm={kernel_norm:.3f}")

        finally:
            # Always cleanup buildStamps context after fitKernel,
            # since fitKernel may call check_again → fillStamp which needs mRData.
            lib = _hotpants_ffi.get_library()
            lib.cleanupBuildStampsContext()

        # Free stamp memory
        n_total_stamps = (
            layout.n_regions_x
            * layout.n_regions_y
            * layout.stamps_per_region_x
            * layout.stamps_per_region_y
        )
        _core.free_stamps(stamps, n_total_stamps)

        return KernelSolution(
            chi2=chi2,
            kernel_norm=kernel_norm,
            mean_sigma=mean_sigma,
            scatter_sigma=scatter_sigma,
            n_skipped_stamps=n_skipped,
            kernel_coefficients=kernel_coeffs,
        )


def spatial_convolve(
    image: np.ndarray,
    kernel_solution: KernelSolution,
    config: KernelConfig | None = None,
    output: np.ndarray | None = None,
    verbose: int = 0,
) -> np.ndarray:
    """
    Apply spatially-varying kernel as convolution to image.

    Evaluates the fitted kernel at every pixel and convolves the image:
        output(x,y) = image(x,y) ⊗ K(x,y)

    where K(x,y) is reconstructed from the kernel coefficients:
        K(x,y) = Σᵢ cᵢ(x,y) * φᵢ(x,y)

    Args:
        image: Input image (numpy array, float32, 2D).

        kernel_solution: KernelSolution from fit_kernel().

        config: KernelConfig instance (default: KernelConfig()).
            Must match the config used in fit_kernel().

        output: Pre-allocated output array (float32, same shape as image).
            If None, a new array is allocated. If provided, must be
            C-contiguous and same shape/dtype as image.

        verbose: Verbosity level (0=silent, 1=progress, 2=debug).

    Returns:
        Output array (float32, same shape as input).
        If output was provided, returns the same array (filled in-place).

    Raises:
        TypeError: if image is not float32
        ValueError: if image shape doesn't match kernel_solution

    Example:
        difference = spatial_convolve(science, kernel_solution)
    """
    if config is None:
        config = KernelConfig()

    validate_image_array(image, "image")

    if output is None:
        output = allocate_array(image.shape, dtype=np.float32)
    else:
        validate_image_array(output, "output")
        if output.shape != image.shape:
            msg = f"Output shape {output.shape} doesn't match input {image.shape}"
            raise ValueError(msg)

    ny, nx = image.shape

    # Make C-contiguous if needed
    if not image.flags["C_CONTIGUOUS"]:
        image = np.ascontiguousarray(image)

    config_dict = {
        "hwKernel": config.kernel_half_width,
        "kerOrder": config.kernel_order,
        "bgOrder": config.bg_order,
        "verbose": verbose,
    }

    with global_state(config_dict):
        logger.debug("Applying spatially-varying convolution...")
        output = _core.spatial_convolve_c(image, kernel_solution.kernel_coefficients, output)

        # Add the fitted background polynomial, matching CLI behavior:
        #   oRData += get_background(k, l, kernelSol)  (main.c lines 1201-1203)
        # The kernel fit models: science ≈ template⊗K + background_poly(x,y).
        # Callers compute D = science - output, so output must include the background
        # term so that D matches the CLI difference image.
        #
        # C get_background() uses: backgroundComponentOffset = (nCompKer-1)*nComp + 1
        # and solutionIdx starting at 1, so first background coefficient is at
        # kernelSol[(nCompKer-1)*nComp + 2].
        n_comp_ker = _hotpants_ffi.get_global_int("nCompKer")
        n_comp = _hotpants_ffi.get_global_int("nComp")
        if n_comp_ker > 0 and n_comp > 0:
            bg_start = (n_comp_ker - 1) * n_comp + 2
            n_bg_terms = (config.bg_order + 1) * (config.bg_order + 2) // 2
            if bg_start + n_bg_terms <= len(kernel_solution.kernel_coefficients):
                bg_coeffs = kernel_solution.kernel_coefficients[bg_start:bg_start + n_bg_terms]
                # Normalise pixel coordinates to [-1, 1] — matches C get_background()
                norm_x = (np.arange(nx, dtype=np.float64) - 0.5 * nx) / (0.5 * nx)
                norm_y = (np.arange(ny, dtype=np.float64) - 0.5 * ny) / (0.5 * ny)
                bg = np.zeros((ny, nx), dtype=np.float64)
                coeff_idx = 0
                for deg_x in range(config.bg_order + 1):
                    for deg_y in range(config.bg_order - deg_x + 1):
                        bg += bg_coeffs[coeff_idx] * np.outer(norm_y ** deg_y, norm_x ** deg_x)
                        coeff_idx += 1
                output += bg.astype(np.float32)

    logger.debug(f"Convolution complete: output shape {output.shape}")

    return output
