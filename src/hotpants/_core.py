"""
Low-level ctypes wrappers and numpy-to-C conversions.

This module provides:
- Array conversion utilities (numpy ↔ C pointers)
- Global state management (set/restore configuration variables)
- Low-level wrapper functions around core C functions
- Error handling and validation

Functions here are not intended for direct use by Python users;
they are internal utilities called by the high-level API in api.py.
"""

import ctypes
from collections.abc import Generator
from contextlib import contextmanager

import numpy as np
from numpy.typing import DTypeLike

from . import _hotpants_ffi


def validate_image_array(arr: np.ndarray, name: str = "image") -> None:
    """
    Validate image array for use with C code.

    Checks:
    - dtype is float32
    - array is 2D
    - no NaN or Inf values

    Note: Non-C-contiguous arrays are acceptable and will be converted
    automatically during processing.

    Args:
        arr: numpy array to validate
        name: name for error messages

    Raises:
        TypeError: if dtype is not float32
        ValueError: if array shape or values are invalid
    """
    if arr.dtype != np.float32:
        msg = f"{name}: expected float32, got {arr.dtype}. Convert with arr.astype(np.float32)."
        raise TypeError(msg)
    if arr.ndim != 2:  # noqa: PLR2004
        msg_0 = f"{name}: expected 2D array, got shape {arr.shape}"
        raise ValueError(msg_0)
    if not np.isfinite(arr).all():
        msg_1 = f"{name}: contains NaN or Inf values. Check for bad pixels or invalid thresholds."
        raise ValueError(msg_1)


def array_to_cptr(arr: np.ndarray):  # noqa: ANN202
    """
    Convert numpy array to C float pointer (borrowed reference).

    Creates a C float* that points to the numpy array's data.
    The numpy array must stay in scope during the C function call.

    Args:
        arr: numpy array (float32, 2D, C-contiguous)

    Returns:
        ctypes pointer to array data
    """
    if arr.dtype != np.float32:
        msg = f"Expected float32, got {arr.dtype}"
        raise TypeError(msg)
    if not arr.flags["C_CONTIGUOUS"]:
        arr = np.ascontiguousarray(arr)

    # Cast numpy array pointer to ctypes void*
    # This is a borrowed reference; numpy array owns the memory
    return ctypes.cast(arr.ctypes.data, ctypes.c_void_p)


def array_to_double_cptr(arr: np.ndarray) -> ctypes.c_void_p:
    """
    Convert numpy array (float32 or float64) to C double pointer.

    If input is float32, creates a temporary float64 copy.
    The resulting pointer is valid only during C function execution.

    Args:
        arr: numpy array (float32 or float64, 2D, C-contiguous)

    Returns:
        ctypes pointer to double data
    """
    if arr.dtype == np.float32:
        arr = arr.astype(np.float64)
    elif arr.dtype != np.float64:
        msg = f"Expected float32/float64, got {arr.dtype}"
        raise TypeError(msg)

    if not arr.flags["C_CONTIGUOUS"]:
        arr = np.ascontiguousarray(arr)

    return ctypes.cast(arr.ctypes.data, ctypes.c_void_p)


def allocate_array(shape: tuple[int, int], dtype: DTypeLike = np.float32) -> np.ndarray:
    """
    Allocate a numpy array suitable for C code.

    Arrays are always C-contiguous (row-major), required by C code.

    Args:
        shape: tuple (ny, nx)
        dtype: data type (default: float32)

    Returns:
        numpy array, C-contiguous, zero-initialized
    """
    return np.zeros(shape, dtype=dtype, order="C")


@contextmanager
def global_state(config_dict: dict[str, int | float]) -> Generator:
    """
    Context manager to temporarily set global configuration variables.

    Sets globals before C function calls, restores on exit.

    Args:
        config_dict: dict mapping global variable names to values
                     e.g., {'hwKernel': 10, 'kerOrder': 2}

    Yields:
        None

    Example:
        with global_state({'hwKernel': 10, 'kerOrder': 2}):
            fitKernel(...)  # Uses configured globals
    """
    # Identify int vs float variables
    int_vars = {
        "hwKernel",
        "kerOrder",
        "bgOrder",
        "nKSStamps",
        "hwKSStamp",
        "nRegX",
        "nRegY",
        "nStampX",
        "nStampY",
        "useFullSS",
        "verbose",
        "nThread",
        "nCompKer",
        "nComp",
        "nCompBG",
    }
    float_vars = {
        "kerFitThresh",
        "scaleFitThresh",
        "tUThresh",
        "tLThresh",
        "iUThresh",
        "iLThresh",
        "tGain",
        "iGain",
        "tRdnoise",
        "iRdnoise",
        "tPedestal",
        "iPedestal",
    }

    # Save original values
    saved = {}
    for key in config_dict:
        if key in int_vars:
            saved[key] = _hotpants_ffi.get_global_int(key)
        elif key in float_vars:
            saved[key] = _hotpants_ffi.get_global_float(key)
        else:
            msg = f"Unknown global variable: {key}"
            raise AttributeError(msg)

    # Set new values
    try:
        for key, val in config_dict.items():
            if key in int_vars:
                _hotpants_ffi.set_global_int(key, int(val))
            elif key in float_vars:
                _hotpants_ffi.set_global_float(key, float(val))
        yield
    finally:
        # Restore original values
        for key, val in saved.items():
            if key in int_vars:
                _hotpants_ffi.set_global_int(key, int(val))
            elif key in float_vars:
                _hotpants_ffi.set_global_float(key, float(val))


# =====================================================================
# Stamp Structure Definition (ctypes)
# =====================================================================


class StampStruct(ctypes.Structure):
    """
    ctypes definition of stamp_struct from hotpants_api.h.

    Represents a single cutout region (stamp) where the kernel is fit.
    Fields map directly to the C structure used by buildStamps() and fitKernel().

    Reference: Alard & Lupton (1998), Section 2.1
    """

    # Define fields before the class is fully constructed (for self-referential pointers)
    pass


# Populate fields after declaration (to handle double** pointers)
StampStruct._fields_ = [
    ("x0", ctypes.c_int),
    ("y0", ctypes.c_int),
    ("x", ctypes.c_int),
    ("y", ctypes.c_int),
    ("nx", ctypes.c_int),
    ("ny", ctypes.c_int),
    ("xss", ctypes.POINTER(ctypes.c_int)),
    ("yss", ctypes.POINTER(ctypes.c_int)),
    ("nss", ctypes.c_int),
    ("sscnt", ctypes.c_int),
    ("vectors", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ("krefArea", ctypes.POINTER(ctypes.c_double)),
    ("mat", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ("scprod", ctypes.POINTER(ctypes.c_double)),
    ("sum", ctypes.c_double),
    ("mean", ctypes.c_double),
    ("median", ctypes.c_double),
    ("mode", ctypes.c_double),
    ("sd", ctypes.c_double),
    ("fwhm", ctypes.c_double),
    ("lfwhm", ctypes.c_double),
    ("chi2", ctypes.c_double),
    ("norm", ctypes.c_double),
    ("diff", ctypes.c_double),
]


# =====================================================================
# Stamp Management Wrappers
# =====================================================================


def allocate_stamps(n_stamps: int) -> ctypes.POINTER(StampStruct):
    """
    Allocate stamp array via C library.

    Args:
        n_stamps: number of stamps to allocate

    Returns:
        ctypes pointer to stamp array

    Raises:
        RuntimeError: if allocation fails in C code
    """
    lib = _hotpants_ffi.get_library()
    allocate_c = lib.allocateStamps
    allocate_c.argtypes = [ctypes.POINTER(StampStruct), ctypes.c_int]
    allocate_c.restype = ctypes.c_int

    # Allocate Python array
    stamps = (StampStruct * n_stamps)()

    # Call C function to initialize internal fields
    result = allocate_c(stamps, n_stamps)
    if result != 0:
        msg = f"C allocateStamps failed with return code {result}"
        raise RuntimeError(msg)

    return stamps


def free_stamps(stamps: ctypes.POINTER(StampStruct), n_stamps: int) -> None:
    """
    Free stamp array memory allocated by C library.

    Args:
        stamps: pointer to stamp array
        n_stamps: number of stamps
    """
    lib = _hotpants_ffi.get_library()
    free_c = lib.freeStampMem
    free_c.argtypes = [ctypes.POINTER(StampStruct), ctypes.c_int]
    free_c.restype = None

    free_c(stamps, n_stamps)


# =====================================================================
# Image Statistics Wrapper
# =====================================================================


def get_stamp_stats(
    data_array: np.ndarray, verbose: int = 0, n_thread: int = 1, n_comp: int = 3
) -> dict[str, float | int]:
    """
    Compute robust image statistics via histogram analysis.

    Estimates mean, median, mode, standard deviation, and FWHM
    using sigma-clipping on a histogram.

    Args:
        data_array: input pixel array (float32, 2D)
        verbose: verbosity level (0=silent, 1=progress, 2=debug)
        n_thread: number of threads
        n_comp: number of kernel components (for internal use)

    Returns:
        dict with keys: 'mean', 'median', 'mode', 'sd', 'fwhm', 'lfwhm'

    Raises:
        ValueError: if data is invalid
    """
    validate_image_array(data_array, "data")

    ny, nx = data_array.shape
    data_ptr = array_to_cptr(data_array)

    # Allocate output variables as ctypes
    mean_val = ctypes.c_double()
    median_val = ctypes.c_double()
    mode_val = ctypes.c_double()
    sd_val = ctypes.c_double()
    fwhm_val = ctypes.c_double()
    lfwhm_val = ctypes.c_double()

    # Get the C function
    lib = _hotpants_ffi.get_library()
    get_stamp_stats3 = lib.getStampStats3
    get_stamp_stats3.argtypes = [
        ctypes.c_void_p,  # float *data
        ctypes.c_int,  # nx
        ctypes.c_int,  # ny
        ctypes.c_int,  # nsy
        ctypes.c_int,  # stat_type
        ctypes.POINTER(ctypes.c_double),  # *mean
        ctypes.POINTER(ctypes.c_double),  # *median
        ctypes.POINTER(ctypes.c_double),  # *mode
        ctypes.POINTER(ctypes.c_double),  # *sd
        ctypes.POINTER(ctypes.c_double),  # *fwhm
        ctypes.POINTER(ctypes.c_double),  # *lfwhm
        ctypes.c_int,  # verbose
        ctypes.c_int,  # nThread
        ctypes.c_int,  # nComp
    ]
    get_stamp_stats3.restype = ctypes.c_int

    # Call C function
    n_used = get_stamp_stats3(
        data_ptr,
        nx,
        ny,
        ny,  # nx, ny, nsy (row count)
        0,  # stat_type (unused)
        ctypes.byref(mean_val),
        ctypes.byref(median_val),
        ctypes.byref(mode_val),
        ctypes.byref(sd_val),
        ctypes.byref(fwhm_val),
        ctypes.byref(lfwhm_val),
        verbose,
        n_thread,
        n_comp,
    )

    return {
        "n_used": n_used,
        "mean": mean_val.value,
        "median": median_val.value,
        "mode": mode_val.value,
        "sd": sd_val.value,
        "fwhm": fwhm_val.value,
        "lfwhm": lfwhm_val.value,
    }


# =====================================================================
# Polynomial Basis Term Count
# =====================================================================


def num_polynomial_terms(order: int) -> int:
    """
    Compute number of polynomial basis terms up to given order.

    For 2D polynomials up to degree 'order':
    N = (order+1)*(order+2)/2

    This counts all monomials x^i * y^j where i,j >= 0 and i+j <= order.

    Example:
        order=2 → 6 terms: {1, x, y, x², xy, y²}

    Reference: Alard & Lupton (1998), Eq. 3

    Args:
        order: polynomial degree (0, 1, 2, ...)

    Returns:
        int: number of terms
    """
    return (order + 1) * (order + 2) // 2


def get_kernel_info() -> dict[str, int]:
    """
    Query current kernel configuration from globals.

    Returns:
        dict with keys:
        - 'kernel_components': nCompKer (number of Gaussian basis functions)
        - 'kernel_order': kerOrder (spatial polynomial order)
        - 'bg_order': bgOrder (background polynomial order)
        - 'kernel_terms': nComp (spatial polynomial terms)
        - 'bg_terms': nCompBG (background polynomial terms)
    """
    return {
        "kernel_components": _hotpants_ffi.get_global_int("nCompKer"),
        "kernel_order": _hotpants_ffi.get_global_int("kerOrder"),
        "bg_order": _hotpants_ffi.get_global_int("bgOrder"),
        "kernel_terms": _hotpants_ffi.get_global_int("nComp"),
        "bg_terms": _hotpants_ffi.get_global_int("nCompBG"),
    }


# =====================================================================
# Core Algorithm Wrappers
# =====================================================================


def build_stamps(
    template: np.ndarray,
    science: np.ndarray,
    n_regions_x: int = 1,
    n_regions_y: int = 1,
    stamps_per_region_x: int = 10,
    stamps_per_region_y: int = 10,
) -> ctypes.POINTER(StampStruct):
    """
    Build stamp grid and compute statistics.

    Divides image into regions and stamp grids, identifies PSF centers,
    and computes stamp statistics via C buildStamps wrapper functions.

    Algorithm:
    1. Initialize buildStamps context with full images
    2. For each region: call C wrapper buildStampsRegion
    3. Clean up context at the end

    Args:
        template: template image (float32, shape ny×nx)
        science: science image (float32, same shape as template)
        n_regions_x, n_regions_y: number of regions per axis
        stamps_per_region_x, stamps_per_region_y: stamps per region per axis

    Returns:
        Pointer to allocated stamp array (total: nRegX×nRegY×nStampX×nStampY)

    Reference: Alard & Lupton (1998) Section 2.1; src/main.c lines 746-985
    """
    validate_image_array(template, "template")
    validate_image_array(science, "science")

    if template.shape != science.shape:
        msg = f"Template and science images must have same shape, got {template.shape} vs {science.shape}"
        raise ValueError(msg)

    ny, nx = science.shape

    # Ensure images are C-contiguous
    if not template.flags["C_CONTIGUOUS"]:
        template = np.ascontiguousarray(template)
    if not science.flags["C_CONTIGUOUS"]:
        science = np.ascontiguousarray(science)

    # Total stamps across all regions
    n_stamps = n_regions_x * n_regions_y * stamps_per_region_x * stamps_per_region_y

    # Allocate primary stamp array
    stamps = allocate_stamps(n_stamps)
    if not stamps:
        msg = "Failed to allocate stamp array"
        raise RuntimeError(msg)

    try:
        # Get wrapper functions from the C library
        lib = _hotpants_ffi.get_library()
        init_context = lib.initBuildStampsContext
        build_region = lib.buildStampsRegion
        cleanup_context = lib.cleanupBuildStampsContext

        # Set argtypes for wrapper functions
        init_context.argtypes = [
            ctypes.c_void_p,  # float *template
            ctypes.c_void_p,  # float *science
            ctypes.c_int,  # ny
            ctypes.c_int,  # nx
            ctypes.c_int,  # n_regions_x
            ctypes.c_int,  # n_regions_y
            ctypes.c_int,  # stamps_per_region_x
            ctypes.c_int,  # stamps_per_region_y
        ]
        init_context.restype = ctypes.c_int

        build_region.argtypes = [
            ctypes.c_int,  # region_x
            ctypes.c_int,  # region_y
            ctypes.POINTER(ctypes.c_void_p),  # *out_template_region
            ctypes.POINTER(ctypes.c_void_p),  # *out_science_region
        ]
        build_region.restype = ctypes.c_int

        # Get the C buildStamps function for direct calls
        build_c = lib.buildStamps
        build_c.argtypes = [
            ctypes.c_int,  # sXMin
            ctypes.c_int,  # sXMax
            ctypes.c_int,  # sYMin
            ctypes.c_int,  # sYMax
            ctypes.POINTER(ctypes.c_int),  # *niS (science stamp counter)
            ctypes.POINTER(ctypes.c_int),  # *ntS (template stamp counter)
            ctypes.c_int,  # getCenters (0=use center, 1=find bright centers)
            ctypes.c_int,  # rXBMin (region x boundary minimum)
            ctypes.c_int,  # rYBMin (region y boundary minimum)
            ctypes.POINTER(StampStruct),  # *ciStamps (science stamps)
            ctypes.POINTER(StampStruct),  # *ctStamps (template stamps)
            ctypes.c_void_p,  # float *iRData (science region data)
            ctypes.c_void_p,  # float *tRData (template region data)
            ctypes.c_float,  # hardX (hard-coded x center, 0.0 to skip)
            ctypes.c_float,  # hardY (hard-coded y center, 0.0 to skip)
        ]
        build_c.restype = None

        cleanup_context.argtypes = []
        cleanup_context.restype = None

        # Get C function pointers for images
        template_ptr = array_to_cptr(template)
        science_ptr = array_to_cptr(science)

        # Set required C globals for buildStamps
        _hotpants_ffi.set_global_int("nKSStamps", 3)  # number of kernel test substamps
        _hotpants_ffi.set_global_int("hwKSStamp", 15)  # half size of kernel substamp
        _hotpants_ffi.set_global_float("tUKThresh", _hotpants_ffi.get_global_float("tUThresh"))

        # Get kernel parameters
        hw_kernel = _hotpants_ffi.get_global_int("hwKernel")
        fw_kernel = hw_kernel * 2 + 1
        _hotpants_ffi.set_global_int("fwKernel", fw_kernel)
        _hotpants_ffi.set_global_int("fwStamp", fw_kernel + 10)

        # Initialize buildStamps context with full images
        ret = init_context(
            template_ptr,  # template image
            science_ptr,  # science image
            ny,  # full image height
            nx,  # full image width
            n_regions_x,  # number of regions x
            n_regions_y,  # number of regions y
            stamps_per_region_x,  # stamps per region x
            stamps_per_region_y,  # stamps per region y
        )

        if ret != 0:
            msg = "Failed to initialize buildStamps context"
            raise RuntimeError(msg)

        try:
            # Process each region
            stamps_per_region = stamps_per_region_x * stamps_per_region_y
            global_stamp_idx = 0

            # Calculate region dimensions
            region_width = nx // n_regions_x
            region_height = ny // n_regions_y

            for region_y in range(n_regions_y):
                for region_x in range(n_regions_x):
                    # Get region boundaries in full image coordinates
                    r_x_min = region_x * region_width
                    r_x_max = min((region_x + 1) * region_width - 1, nx - 1)
                    r_y_min = region_y * region_height
                    r_y_max = min((region_y + 1) * region_height - 1, ny - 1)

                    # Add kernel border for edge handling
                    r_x_b_min = max(0, r_x_min - hw_kernel)
                    r_x_b_max = min(nx - 1, r_x_max + hw_kernel)
                    r_y_b_min = max(0, r_y_min - hw_kernel)
                    r_y_b_max = min(ny - 1, r_y_max + hw_kernel)

                    # Call wrapper to extract region data and setup globals
                    region_t_ptr = ctypes.c_void_p()
                    region_i_ptr = ctypes.c_void_p()
                    ret = build_region(
                        region_x,  # region x coordinate
                        region_y,  # region y coordinate
                        ctypes.byref(region_t_ptr),  # output template region pointer
                        ctypes.byref(region_i_ptr),  # output science region pointer
                    )

                    if ret != 0:
                        msg = f"Failed to extract region ({region_x}, {region_y})"
                        raise RuntimeError(msg)

                    # Now call buildStamps directly with extracted region data
                    # Iterate over stamps within this region
                    for stamp_y in range(stamps_per_region_y):
                        for stamp_x in range(stamps_per_region_x):
                            # Calculate stamp bounds relative to region borders
                            stamp_width = (r_x_b_max - r_x_b_min + 1) // stamps_per_region_x
                            stamp_height = (r_y_b_max - r_y_b_min + 1) // stamps_per_region_y

                            s_x_min = stamp_x * stamp_width
                            s_y_min = stamp_y * stamp_height
                            s_x_max = min(s_x_min + stamp_width - 1, (r_x_b_max - r_x_b_min))
                            s_y_max = min(s_y_min + stamp_height - 1, (r_y_b_max - r_y_b_min))

                            # Stamp counters
                            ni_s = ctypes.c_int(0)
                            nt_s = ctypes.c_int(0)

                            # Call buildStamps with region data
                            build_c(
                                s_x_min,  # sXMin
                                s_x_max,  # sXMax
                                s_y_min,  # sYMin
                                s_y_max,  # sYMax
                                ctypes.byref(ni_s),  # *niS
                                ctypes.byref(nt_s),  # *ntS
                                1,  # getCenters: automatically find bright centers
                                r_x_b_min,  # rXBMin: region origin in full image
                                r_y_b_min,  # rYBMin: region origin in full image
                                ctypes.byref(stamps[global_stamp_idx]),  # *ciStamps
                                ctypes.byref(stamps[global_stamp_idx]),  # *ctStamps
                                region_i_ptr,  # float *iRData (region-sized)
                                region_t_ptr,  # float *tRData (region-sized)
                                0.0,  # hardX
                                0.0,  # hardY
                            )

                            global_stamp_idx += 1

        finally:
            # Always clean up context
            cleanup_context()

    except Exception:
        free_stamps(stamps, n_stamps)
        raise

    return stamps


def fit_kernel_c(
    stamps: ctypes.POINTER(StampStruct),
    template: np.ndarray,
    science: np.ndarray,
    noise: np.ndarray | None = None,
    n_kernel_components: int = 3,
    n_spatial_terms: int = 6,
) -> tuple[float, float, float, float, int, np.ndarray]:
    """
    Fit kernel via least-squares.

    Args:
        stamps: stamp array pointer
        template: template image
        science: science image
        noise: noise image (optional)
        n_kernel_components: number of kernel basis components
        n_spatial_terms: number of spatial polynomial terms

    Returns:
        Tuple: (chi2, kernel_norm, mean_sigma, scatter_sigma, n_skipped, coefficients)
    """
    lib = _hotpants_ffi.get_library()
    fit_c = lib.fitKernel
    fit_c.argtypes = [
        ctypes.c_void_p,  # stamp_struct *stamps
        ctypes.c_void_p,  # float *imRef
        ctypes.c_void_p,  # float *imConv
        ctypes.c_void_p,  # float *imNoise
        ctypes.c_void_p,  # double *kernel_coeffs
        ctypes.POINTER(ctypes.c_double),  # *meansig
        ctypes.POINTER(ctypes.c_double),  # *scatter
        ctypes.POINTER(ctypes.c_int),  # *n_skipped
    ]
    fit_c.restype = None

    # Allocate output variables
    n_coeffs = n_kernel_components * n_spatial_terms
    kernel_coeffs = np.zeros(n_coeffs, dtype=np.float64)
    meansig = ctypes.c_double()
    scatter = ctypes.c_double()
    n_skipped = ctypes.c_int()

    # Convert to C pointers
    tpl_ptr = array_to_cptr(template)
    sci_ptr = array_to_cptr(science)
    noise_ptr = array_to_cptr(noise) if noise is not None else None
    coeffs_ptr = array_to_double_cptr(kernel_coeffs)

    # Call fitKernel
    fit_c(
        stamps,
        tpl_ptr,
        sci_ptr,
        noise_ptr,
        coeffs_ptr,
        ctypes.byref(meansig),
        ctypes.byref(scatter),
        ctypes.byref(n_skipped),
    )

    # TODO: fitKernel in C doesn't return chi2 directly; need to compute or extract from stamps
    chi2 = 0.0  # placeholder

    return (chi2, 1.0, meansig.value, scatter.value, n_skipped.value, kernel_coeffs)


def spatial_convolve_c(
    image: np.ndarray, kernel_coeffs: np.ndarray, output: np.ndarray | None = None
) -> np.ndarray:
    """
    Apply spatially-varying kernel convolution.

    Args:
        image: input image
        kernel_coeffs: kernel coefficients from fit_kernel_c
        output: pre-allocated output array (optional)

    Returns:
        Convolved image
    """
    if output is None:
        output = allocate_array(image.shape, dtype=np.float32)

    lib = _hotpants_ffi.get_library()
    convolve_c = lib.spatial_convolve
    convolve_c.argtypes = [
        ctypes.c_void_p,  # float *image
        ctypes.c_void_p,  # float **var_image
        ctypes.c_int,  # ny
        ctypes.c_int,  # nx
        ctypes.c_void_p,  # double *kernel_coeffs
        ctypes.c_void_p,  # float *output
        ctypes.POINTER(ctypes.c_int),  # *conv_method
    ]
    convolve_c.restype = None

    ny, nx = image.shape

    # Convert to C pointers
    img_ptr = array_to_cptr(image)
    out_ptr = array_to_cptr(output)
    coeffs_ptr = array_to_double_cptr(kernel_coeffs)

    # Get convolution method (0=direct, 1=FFT)
    conv_method = ctypes.c_int(1)  # Use FFT if available

    # Call spatial_convolve
    convolve_c(img_ptr, None, ny, nx, coeffs_ptr, out_ptr, ctypes.byref(conv_method))

    return output
