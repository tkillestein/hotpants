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

from dataclasses import dataclass
from typing import Optional
import numpy as np

from . import _core
from ._core import global_state, validate_image_array, allocate_array


# =====================================================================
# Configuration Dataclasses
# =====================================================================

@dataclass
class KernelConfig:
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
    kernel_half_width: int = 15
    kernel_order: int = 2
    bg_order: int = 1
    fit_threshold: float = 20.0
    scale_fit_threshold: float = 0.5
    n_ks_stamps: int = 3
    hw_ks_stamp: int = 10

    def __post_init__(self):
        """Validate configuration parameters."""
        if self.kernel_half_width <= 0:
            raise ValueError("kernel_half_width must be positive")
        if self.kernel_order < 0:
            raise ValueError("kernel_order must be >= 0")
        if self.bg_order < 0:
            raise ValueError("bg_order must be >= 0")
        if self.fit_threshold <= 0:
            raise ValueError("fit_threshold must be positive")
        if not (0 < self.scale_fit_threshold <= 1):
            raise ValueError("scale_fit_threshold must be in (0, 1]")
        if self.n_ks_stamps <= 0:
            raise ValueError("n_ks_stamps must be positive")
        if self.hw_ks_stamp <= 0:
            raise ValueError("hw_ks_stamp must be positive")
        if self.hw_ks_stamp > self.kernel_half_width:
            raise ValueError("hw_ks_stamp must be <= kernel_half_width")


@dataclass
class RegionLayout:
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
    n_regions_x: int = 1
    n_regions_y: int = 1
    stamps_per_region_x: int = 10
    stamps_per_region_y: int = 10
    use_full_substamps: bool = False

    def __post_init__(self):
        """Validate layout parameters."""
        if self.n_regions_x < 1 or self.n_regions_y < 1:
            raise ValueError("n_regions_* must be >= 1")
        if self.stamps_per_region_x < 1 or self.stamps_per_region_y < 1:
            raise ValueError("stamps_per_region_* must be >= 1")


@dataclass
class NoiseThresholds:
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
    template_gain: float = 1.0
    science_gain: float = 1.0
    template_readnoise: float = 0.0
    science_readnoise: float = 0.0
    template_pedestal: float = 0.0
    science_pedestal: float = 0.0

    def __post_init__(self):
        """Validate threshold parameters."""
        if self.template_gain <= 0 or self.science_gain <= 0:
            raise ValueError("gain values must be positive")
        if self.template_readnoise < 0 or self.science_readnoise < 0:
            raise ValueError("readnoise values must be >= 0")
        if self.template_upper_threshold <= self.template_lower_threshold:
            raise ValueError("upper_threshold must be > lower_threshold")


@dataclass
class KernelSolution:
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
    chi2: float
    kernel_norm: float
    mean_sigma: float
    scatter_sigma: float
    n_skipped_stamps: int
    kernel_coefficients: np.ndarray

    def __repr__(self):
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
    noise: Optional[np.ndarray] = None,
    config: Optional[KernelConfig] = None,
    layout: Optional[RegionLayout] = None,
    thresholds: Optional[NoiseThresholds] = None,
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
        raise ValueError(
            f"Image shape mismatch: template {template.shape} != science {science.shape}"
        )

    if noise is not None:
        validate_image_array(noise, "noise")
        if noise.shape != template.shape:
            raise ValueError(
                f"Noise shape mismatch: {noise.shape} != {template.shape}"
            )

    ny, nx = template.shape

    # Convert to C-contiguous if needed
    if not template.flags['C_CONTIGUOUS']:
        template = np.ascontiguousarray(template)
    if not science.flags['C_CONTIGUOUS']:
        science = np.ascontiguousarray(science)
    if noise is not None and not noise.flags['C_CONTIGUOUS']:
        noise = np.ascontiguousarray(noise)

    # Set global configuration
    config_dict = {
        'hwKernel': config.kernel_half_width,
        'kerOrder': config.kernel_order,
        'bgOrder': config.bg_order,
        'kerFitThresh': config.fit_threshold,
        'scaleFitThresh': config.scale_fit_threshold,
        'nKSStamps': config.n_ks_stamps,
        'hwKSStamp': config.hw_ks_stamp,
        'nRegX': layout.n_regions_x,
        'nRegY': layout.n_regions_y,
        'nStampX': layout.stamps_per_region_x,
        'nStampY': layout.stamps_per_region_y,
        'useFullSS': 1 if layout.use_full_substamps else 0,
        'tUThresh': thresholds.template_upper_threshold,
        'tLThresh': thresholds.template_lower_threshold,
        'iUThresh': thresholds.science_upper_threshold,
        'iLThresh': thresholds.science_lower_threshold,
        'tGain': thresholds.template_gain,
        'iGain': thresholds.science_gain,
        'tRdnoise': thresholds.template_readnoise,
        'iRdnoise': thresholds.science_readnoise,
        'tPedestal': thresholds.template_pedestal,
        'iPedestal': thresholds.science_pedestal,
        'verbose': verbose,
        'nThread': n_thread,
    }

    with global_state(config_dict):
        # Placeholder: actual implementation would call buildStamps, fitKernel, etc.
        # For now, return a dummy result to verify the infrastructure works
        kernel_info = _core.get_kernel_info()
        n_terms = _core.num_polynomial_terms(config.kernel_order)
        n_kernel_comps = kernel_info['kernel_components']

        # Create dummy kernel coefficients
        kernel_coeffs = np.ones(n_kernel_comps * n_terms, dtype=np.float64)

        return KernelSolution(
            chi2=0.0,
            kernel_norm=1.0,
            mean_sigma=10.0,
            scatter_sigma=1.0,
            n_skipped_stamps=0,
            kernel_coefficients=kernel_coeffs,
        )


def spatial_convolve(
    image: np.ndarray,
    kernel_solution: KernelSolution,
    config: Optional[KernelConfig] = None,
    output: Optional[np.ndarray] = None,
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
            raise ValueError(
                f"Output shape {output.shape} doesn't match input {image.shape}"
            )

    # Placeholder: actual implementation would call spatial_convolve C function
    # For now, copy input to output to verify infrastructure works
    output[:] = image

    return output
