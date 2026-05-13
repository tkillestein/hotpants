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
from loguru import logger
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
        "nCompTotal",
        "nS",
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
# Context Structure Definitions (ctypes)
# =====================================================================


class ConfigStruct(ctypes.Structure):
    """ctypes definition of config_struct from hotpants_context.h."""

    _fields_ = [
        ("kernel_half_width", ctypes.c_int),
        ("kernel_order", ctypes.c_int),
        ("bg_order", ctypes.c_int),
        ("n_ks_stamps", ctypes.c_int),
        ("hw_ks_stamp", ctypes.c_int),
        ("template_upper_threshold", ctypes.c_float),
        ("template_lower_threshold", ctypes.c_float),
        ("science_upper_threshold", ctypes.c_float),
        ("science_lower_threshold", ctypes.c_float),
        ("template_gain", ctypes.c_float),
        ("science_gain", ctypes.c_float),
        ("template_readnoise", ctypes.c_float),
        ("science_readnoise", ctypes.c_float),
        ("template_pedestal", ctypes.c_float),
        ("science_pedestal", ctypes.c_float),
        ("fit_threshold", ctypes.c_float),
        ("scale_fit_threshold", ctypes.c_float),
        ("min_frac_good_stamps", ctypes.c_float),
        ("stat_sig", ctypes.c_float),
        ("ker_sig_reject", ctypes.c_float),
        ("ker_frac_mask", ctypes.c_float),
        ("fill_val", ctypes.c_float),
        ("fill_val_noise", ctypes.c_float),
        ("kf_spread_mask1", ctypes.c_float),
        ("kf_spread_mask2", ctypes.c_float),
        ("n_regions_x", ctypes.c_int),
        ("n_regions_y", ctypes.c_int),
        ("stamps_per_region_x", ctypes.c_int),
        ("stamps_per_region_y", ctypes.c_int),
        ("use_full_ss", ctypes.c_int),
        ("verbose", ctypes.c_int),
        ("n_thread", ctypes.c_int),
        ("force_convolve", ctypes.c_char),
    ]


class GaussianBasisStruct(ctypes.Structure):
    """ctypes definition of gaussian_basis_struct from hotpants_context.h."""

    pass


GaussianBasisStruct._fields_ = [
    ("ngauss", ctypes.c_int),
    ("deg_fixe", ctypes.POINTER(ctypes.c_int)),
    ("sigma_gauss", ctypes.POINTER(ctypes.c_float)),
    ("n_comp_ker", ctypes.c_int),
    ("n_comp", ctypes.c_int),
    ("n_c", ctypes.c_int),
    ("n_comp_bg", ctypes.c_int),
    ("n_bg_vectors", ctypes.c_int),
    ("n_comp_total", ctypes.c_int),
    ("kernel_vec", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ("filter_x", ctypes.POINTER(ctypes.c_double)),
    ("filter_y", ctypes.POINTER(ctypes.c_double)),
]


class KernelContextStruct(ctypes.Structure):
    """ctypes definition of kernel_context_struct from hotpants_context.h."""

    _fields_ = [
        ("fw_kernel", ctypes.c_int),
        ("fw_stamp", ctypes.c_int),
        ("fw_ks_stamp", ctypes.c_int),
        ("kc_step", ctypes.c_int),
        ("basis", GaussianBasisStruct),
        ("is_initialized", ctypes.c_int),
    ]


class RegionContextStruct(ctypes.Structure):
    """ctypes definition of region_context_struct from hotpants_context.h."""

    pass


RegionContextStruct._fields_ = [
    ("pix_x", ctypes.c_int),
    ("pix_y", ctypes.c_int),
    ("mask", ctypes.POINTER(ctypes.c_int)),
    ("temp", ctypes.POINTER(ctypes.c_float)),
    ("temp2", ctypes.POINTER(ctypes.c_float)),
    ("wxy", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ("kernel", ctypes.POINTER(ctypes.c_double)),
    ("kernel_coeffs", ctypes.POINTER(ctypes.c_double)),
    ("check_stack", ctypes.POINTER(ctypes.c_double)),
    ("check_mat", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
    ("check_vec", ctypes.POINTER(ctypes.c_double)),
    ("n_stamps_allocated", ctypes.c_int),
    ("is_initialized", ctypes.c_int),
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
# Context Management Wrappers (NEW API)
# =====================================================================


def create_config_struct() -> ConfigStruct:
    """
    Create a configuration struct with default values.

    Returns:
        ConfigStruct instance initialized with defaults
    """
    lib = _hotpants_ffi.get_library()
    default_config_c = lib.hotpants_default_config
    default_config_c.argtypes = [ctypes.POINTER(ConfigStruct)]
    default_config_c.restype = ctypes.c_int

    config = ConfigStruct()
    result = default_config_c(ctypes.byref(config))
    if result != 0:
        msg = "Failed to initialize default config"
        raise RuntimeError(msg)

    return config


def create_kernel_context(
    config: ConfigStruct, image_nx: int, image_ny: int,
    n_regions_x: int, n_regions_y: int,
    stamps_per_region_x: int, stamps_per_region_y: int,
) -> ctypes.POINTER(KernelContextStruct):
    """
    Create and initialize a kernel context from configuration.

    Args:
        config: configuration struct
        image_nx, image_ny: full image dimensions
        n_regions_x, n_regions_y: region grid
        stamps_per_region_x, stamps_per_region_y: stamp grid

    Returns:
        pointer to allocated kernel context

    Raises:
        RuntimeError: if context creation fails
    """
    lib = _hotpants_ffi.get_library()
    create_ctx_c = lib.hotpants_create_kernel_context
    create_ctx_c.argtypes = [
        ctypes.POINTER(ConfigStruct),
        ctypes.c_int, ctypes.c_int,
        ctypes.c_int, ctypes.c_int,
        ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.POINTER(KernelContextStruct)),
    ]
    create_ctx_c.restype = ctypes.c_int

    ctx_ptr = ctypes.POINTER(KernelContextStruct)()
    result = create_ctx_c(
        ctypes.byref(config),
        image_nx, image_ny,
        n_regions_x, n_regions_y,
        stamps_per_region_x, stamps_per_region_y,
        ctypes.byref(ctx_ptr),
    )
    if result != 0:
        msg = "Failed to create kernel context"
        raise RuntimeError(msg)

    return ctx_ptr


def cleanup_kernel_context(ctx: ctypes.POINTER(KernelContextStruct)) -> None:
    """
    Free kernel context memory.

    Args:
        ctx: pointer to kernel context (may be None)
    """
    if not ctx:
        return

    lib = _hotpants_ffi.get_library()
    cleanup_ctx_c = lib.hotpants_cleanup_kernel_context
    cleanup_ctx_c.argtypes = [ctypes.POINTER(KernelContextStruct)]
    cleanup_ctx_c.restype = None

    cleanup_ctx_c(ctx)


def create_region_context(
    kctx: ctypes.POINTER(KernelContextStruct),
    config: ConfigStruct,
    n_stamps: int,
) -> ctypes.POINTER(RegionContextStruct):
    """
    Create and initialize a region context.

    Args:
        kctx: kernel context (must be valid)
        config: configuration struct
        n_stamps: number of stamps for this region

    Returns:
        pointer to allocated region context

    Raises:
        RuntimeError: if context creation fails
    """
    lib = _hotpants_ffi.get_library()
    create_rctx_c = lib.hotpants_create_region_context
    create_rctx_c.argtypes = [
        ctypes.POINTER(KernelContextStruct),
        ctypes.POINTER(ConfigStruct),
        ctypes.c_int,
        ctypes.POINTER(ctypes.POINTER(RegionContextStruct)),
    ]
    create_rctx_c.restype = ctypes.c_int

    rctx_ptr = ctypes.POINTER(RegionContextStruct)()
    result = create_rctx_c(
        kctx,
        ctypes.byref(config),
        n_stamps,
        ctypes.byref(rctx_ptr),
    )
    if result != 0:
        msg = "Failed to create region context"
        raise RuntimeError(msg)

    return rctx_ptr


def cleanup_region_context(ctx: ctypes.POINTER(RegionContextStruct)) -> None:
    """
    Free region context memory.

    Args:
        ctx: pointer to region context (may be None)
    """
    if not ctx:
        return

    lib = _hotpants_ffi.get_library()
    cleanup_rctx_c = lib.hotpants_cleanup_region_context
    cleanup_rctx_c.argtypes = [ctypes.POINTER(RegionContextStruct)]
    cleanup_rctx_c.restype = None

    cleanup_rctx_c(ctx)


# =====================================================================
# Context-Aware Algorithm Wrappers (v2)
# =====================================================================


def fitkernel_v2(
    stamps: ctypes.POINTER(StampStruct),
    n_stamps: int,
    im_ref: np.ndarray,
    im_conv: np.ndarray,
    im_noise: np.ndarray | None,
    config: ConfigStruct,
    kctx: ctypes.POINTER(KernelContextStruct),
    rctx: ctypes.POINTER(RegionContextStruct),
) -> tuple[np.ndarray, float, float, int]:
    """
    Fit kernel using context-aware API (v2).

    Args:
        stamps: pointer to stamp array
        n_stamps: number of stamps
        im_ref: reference image (float32, 2D)
        im_conv: image to convolve (float32, 2D)
        im_noise: noise image or None (float32, 2D)
        config: configuration struct
        kctx: kernel context
        rctx: region context

    Returns:
        tuple of (kernel_coeffs, meansig, scatter, n_skipped)

    Raises:
        RuntimeError: if kernel fitting fails
    """
    lib = _hotpants_ffi.get_library()
    fitkernel_v2_c = lib.fitKernel_v2
    fitkernel_v2_c.argtypes = [
        ctypes.POINTER(StampStruct),
        ctypes.c_int,
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.POINTER(ConfigStruct),
        ctypes.POINTER(KernelContextStruct),
        ctypes.POINTER(RegionContextStruct),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
    ]
    fitkernel_v2_c.restype = ctypes.c_int

    # Convert arrays to C pointers
    validate_image_array(im_ref, "im_ref")
    validate_image_array(im_conv, "im_conv")
    ref_ptr = array_to_cptr(im_ref)
    conv_ptr = array_to_cptr(im_conv)
    noise_ptr = ctypes.c_void_p() if im_noise is None else array_to_cptr(im_noise)

    # Allocate output arrays
    n_comp_total = kctx.contents.basis.n_comp_total
    kernel_coeffs = (ctypes.c_double * (n_comp_total + 1))()
    meansig = ctypes.c_double()
    scatter = ctypes.c_double()
    n_skipped = ctypes.c_int()

    # Call C function
    result = fitkernel_v2_c(
        stamps,
        n_stamps,
        ref_ptr,
        conv_ptr,
        noise_ptr,
        ctypes.byref(config),
        kctx,
        rctx,
        kernel_coeffs,
        ctypes.byref(meansig),
        ctypes.byref(scatter),
        ctypes.byref(n_skipped),
    )

    if result != 0:
        msg = f"fitKernel_v2 failed with return code {result}"
        raise RuntimeError(msg)

    # Convert output to numpy
    coeffs_np = np.array([kernel_coeffs[i] for i in range(n_comp_total + 1)])

    return coeffs_np, meansig.value, scatter.value, n_skipped.value


def spatial_convolve_v2(
    image: np.ndarray,
    kernel_coeffs: np.ndarray,
    config: ConfigStruct,
    kctx: ctypes.POINTER(KernelContextStruct),
    rctx: ctypes.POINTER(RegionContextStruct),
    var_image: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply kernel via convolution using context-aware API (v2).

    Args:
        image: input image (float32, 2D)
        kernel_coeffs: fitted kernel coefficients (float64)
        config: configuration struct
        kctx: kernel context
        rctx: region context
        var_image: variance image or None (float32, 2D)

    Returns:
        tuple of (output_image, conv_mask)

    Raises:
        RuntimeError: if convolution fails
    """
    lib = _hotpants_ffi.get_library()
    convolve_v2_c = lib.spatial_convolve_v2
    convolve_v2_c.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.POINTER(ctypes.c_float)),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ConfigStruct),
        ctypes.POINTER(KernelContextStruct),
        ctypes.POINTER(RegionContextStruct),
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
    ]
    convolve_v2_c.restype = ctypes.c_int

    # Validate inputs
    validate_image_array(image, "image")
    ny, nx = image.shape

    if var_image is not None:
        validate_image_array(var_image, "var_image")

    # Convert arrays
    image_ptr = array_to_cptr(image)
    var_ptr = ctypes.POINTER(ctypes.POINTER(ctypes.c_float))() if var_image is None else None
    coeffs_arr = np.asarray(kernel_coeffs, dtype=np.float64)
    coeffs_ptr = array_to_double_cptr(coeffs_arr)

    # Allocate output
    output = allocate_array((ny, nx), dtype=np.float32)
    output_ptr = array_to_cptr(output)

    conv_mask = (ctypes.c_int * (ny * nx))()

    # Call C function
    result = convolve_v2_c(
        image_ptr,
        var_ptr,
        ny,
        nx,
        ctypes.cast(coeffs_ptr, ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(config),
        kctx,
        rctx,
        output_ptr,
        conv_mask,
    )

    if result != 0:
        msg = f"spatial_convolve_v2 failed with return code {result}"
        raise RuntimeError(msg)

    # Convert mask to numpy
    mask_np = np.array([conv_mask[i] for i in range(ny * nx)], dtype=np.int32)
    mask_np = mask_np.reshape((ny, nx))

    return output, mask_np


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
        "total_components": _hotpants_ffi.get_global_int("nCompTotal"),
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

    try:
        # Get wrapper functions from the C library
        lib = _hotpants_ffi.get_library()
        init_globals = lib.initKernelGlobals
        init_context = lib.initBuildStampsContext
        build_region = lib.buildStampsRegion
        cleanup_context = lib.cleanupBuildStampsContext

        init_globals.argtypes = [
            ctypes.c_int,  # image_nx
            ctypes.c_int,  # image_ny
            ctypes.c_int,  # n_reg_x
            ctypes.c_int,  # n_reg_y
            ctypes.c_int,  # n_stamp_x
            ctypes.c_int,  # n_stamp_y
        ]
        init_globals.restype = ctypes.c_int

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
            ctypes.c_float,  # hardX
            ctypes.c_float,  # hardY
        ]
        build_c.restype = None

        cleanup_context.argtypes = []
        cleanup_context.restype = None

        # Set basis globals before allocating stamps:
        # hwKernel, kerOrder, bgOrder, nKSStamps, hwKSStamp must be in place
        _hotpants_ffi.set_global_int("nKSStamps", _hotpants_ffi.get_global_int("nKSStamps") or 3)
        _hotpants_ffi.set_global_int("hwKSStamp", _hotpants_ffi.get_global_int("hwKSStamp") or 15)
        # Kernel fit upper thresholds: vargs.c sets these from tUThresh/iUThresh;
        # Python callers must mirror that assignment or getPsfCenters treats all
        # pixels as saturated (hiThresh=0 causes dpt >= hiThresh for every pixel).
        _hotpants_ffi.set_global_float("tUKThresh", _hotpants_ffi.get_global_float("tUThresh"))
        _hotpants_ffi.set_global_float("iUKThresh", _hotpants_ffi.get_global_float("iUThresh"))

        # initKernelGlobals computes nCompKer, nC, fwKSStamp, fwKernel, fwStamp
        # (and sets up the Gaussian basis) — must run BEFORE allocateStamps
        ret = init_globals(nx, ny, n_regions_x, n_regions_y,
                           stamps_per_region_x, stamps_per_region_y)
        if ret != 0:
            msg = "Failed to initialize kernel globals"
            raise RuntimeError(msg)

        # Now that nCompKer, nC, fwKSStamp, nKSStamps are set, allocate stamps
        stamps = allocate_stamps(n_stamps)
        if not stamps:
            msg = "Failed to allocate stamp array"
            raise RuntimeError(msg)

        # Get C function pointers for images
        template_ptr = array_to_cptr(template)
        science_ptr = array_to_cptr(science)

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

        # Set up fillStampsForRegion wrapper
        fill_region = lib.fillStampsForRegion
        fill_region.argtypes = [ctypes.c_void_p, ctypes.c_int]
        fill_region.restype = ctypes.c_int

        # Process each region.  For each region:
        #   1. buildStampsRegion extracts the region buffer and sets mRData/rPixX/rPixY
        #   2. buildStamps is called per stamp to find PSF centres (fills xss/yss/nss)
        #   3. fillStampsForRegion fills stamp->vectors/mat/scprod via fillStamp()
        #      while the region context (mRData, buffers) is still active
        # The region context is NOT cleaned up here; the caller (api.py) must
        # call cleanupBuildStampsContext() after fitKernel() has completed.
        stamps_per_region = stamps_per_region_x * stamps_per_region_y
        global_stamp_idx = 0

        # Calculate region dimensions
        region_width = nx // n_regions_x
        region_height = ny // n_regions_y
        hw_kernel = _hotpants_ffi.get_global_int("hwKernel")

        for region_y in range(n_regions_y):
            for region_x in range(n_regions_x):
                region_start_idx = global_stamp_idx

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

                # Now call buildStamps directly with extracted region data.
                # Coordinates mirror main.c lines 931-934: stamp origins are
                # evenly spaced across rPixX/rPixY in full-image coords; stamp
                # extent is fwStamp (not rPixX/nStampX) matching C behaviour.
                r_pix_x = r_x_b_max - r_x_b_min + 1
                r_pix_y = r_y_b_max - r_y_b_min + 1
                fw_stamp = _hotpants_ffi.get_global_int("fwStamp")

                for stamp_y in range(stamps_per_region_y):
                    for stamp_x in range(stamps_per_region_x):
                        # Full-image stamp coordinates (matches main.c)
                        s_x_min = r_x_b_min + stamp_x * r_pix_x // stamps_per_region_x
                        s_y_min = r_y_b_min + stamp_y * r_pix_y // stamps_per_region_y
                        s_x_max = min(s_x_min + fw_stamp - 1, r_x_b_max)
                        s_y_max = min(s_y_min + fw_stamp - 1, r_y_b_max)

                        # Reset stamp state before buildStamps (matches main.c lines 938-943)
                        stamps[global_stamp_idx].nss = 0
                        stamps[global_stamp_idx].sscnt = 0
                        stamps[global_stamp_idx].chi2 = 0.0

                        # Stamp counters (buildStamps uses *niS/*ntS to index ciStamps/ctStamps)
                        ni_s = ctypes.c_int(0)
                        nt_s = ctypes.c_int(0)

                        # Call buildStamps with region data
                        build_c(
                            s_x_min,  # sXMin (full-image coords)
                            s_x_max,  # sXMax (full-image coords)
                            s_y_min,  # sYMin (full-image coords)
                            s_y_max,  # sYMax (full-image coords)
                            ctypes.byref(ni_s),  # *niS
                            ctypes.byref(nt_s),  # *ntS
                            1,  # getCenters: find bright PSF centers
                            r_x_b_min,  # rXBMin: region buffer origin (full-image)
                            r_y_b_min,  # rYBMin: region buffer origin (full-image)
                            ctypes.byref(stamps[global_stamp_idx]),  # *ciStamps
                            ctypes.byref(stamps[global_stamp_idx]),  # *ctStamps
                            region_i_ptr,  # float *iRData (region buffer)
                            region_t_ptr,  # float *tRData (region buffer)
                            0.0,  # hardX
                            0.0,  # hardY
                        )

                        # Mirror main.c: only advance the stamp index when PSF centers
                        # were found (nss > 0). This densely packs valid stamps so
                        # fitKernel processes exactly the same stamps as the CLI.
                        if stamps[global_stamp_idx].nss > 0:
                            global_stamp_idx += 1

                # Fill stamp vectors for this region while mRData + region buffers
                # are still valid (mirrors main.c lines 1099-1119: fillStamp loop).
                region_count = global_stamp_idx - region_start_idx
                region_stamps_ptr = ctypes.addressof(stamps[region_start_idx])
                ret = fill_region(region_stamps_ptr, region_count)
                if ret != 0:
                    msg = f"fillStampsForRegion failed for region ({region_x}, {region_y})"
                    raise RuntimeError(msg)

        # Set nS = number of valid stamps (densely packed at indices 0..nS-1).
        # Mirrors main.c: nS = ntS (only valid template stamps are counted).
        _hotpants_ffi.set_global_int("nS", global_stamp_idx)

        # Region context (mRData, region buffers) is intentionally left active.
        # Caller must call cleanupBuildStampsContext() after fitKernel() completes,
        # because fitKernel() → check_again() → getStampSig/fillStamp access mRData.

    except Exception:
        # On error, clean up context and stamps before re-raising
        try:
            cleanup_context()
        except Exception:
            pass
        try:
            free_stamps(stamps, n_stamps)
        except (UnboundLocalError, Exception):
            pass
        raise

    return stamps


def fit_kernel_c(
    stamps: ctypes.POINTER(StampStruct),
    template: np.ndarray,
    science: np.ndarray,
    noise: np.ndarray | None = None,
    n_coeffs: int = 298,
) -> tuple[float, float, float, float, int, np.ndarray]:
    """
    Fit kernel via least-squares.

    Args:
        stamps: stamp array pointer
        template: template image
        science: science image
        noise: noise image (optional)
        n_coeffs: total coefficient array size = nCompTotal + 1

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

    # Allocate coefficient array: fitKernel writes nCompTotal+1 values
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
    # Get address of stamp array pointer (c_void_p because argtypes specifies c_void_p)
    logger.debug(f"fit_kernel_c: stamps type={type(stamps)}, len={len(stamps)}")

    # Use addressof to get the address of the array
    stamps_addr = ctypes.addressof(stamps)
    stamps_ptr = ctypes.c_void_p(stamps_addr)
    # Note: fitKernel signature is (stamps, imRef, imConv, imNoise, ...)
    # where imRef is the reference/target image and imConv is the image to be convolved.
    # In our case, we fit the kernel to convolve the template to match the science,
    # so: imRef = science (reference), imConv = template (to be convolved)
    fit_c(
        stamps_ptr,
        sci_ptr,  # imRef (science/reference)
        tpl_ptr,  # imConv (template to be convolved)
        noise_ptr,  # imNoise
        coeffs_ptr,
        ctypes.byref(meansig),
        ctypes.byref(scatter),
        ctypes.byref(n_skipped),
    )
    logger.debug("fit_kernel_c: fitKernel returned successfully")

    # Compute actual kernel norm by evaluating the fitted kernel at the image centre.
    # computeKernelNorm sets rPixX/rPixY then calls make_kernel(nx/2, ny/2, coeffs).
    compute_norm_c = lib.computeKernelNorm
    compute_norm_c.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    compute_norm_c.restype = ctypes.c_double
    ny, nx = template.shape
    kernel_norm = compute_norm_c(coeffs_ptr, nx, ny)

    chi2 = 0.0  # fitKernel doesn't return chi2 directly; placeholder

    return (chi2, kernel_norm, meansig.value, scatter.value, n_skipped.value, kernel_coeffs)


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
        ctypes.c_void_p,                    # float *image
        ctypes.POINTER(ctypes.c_void_p),    # float **var_image (pointer to nullable ptr)
        ctypes.c_int,                        # xSize (= nx, image width)
        ctypes.c_int,                        # ySize (= ny, image height)
        ctypes.c_void_p,                    # double *kernel_coeffs
        ctypes.c_void_p,                    # float *output (cRdata)
        ctypes.c_void_p,                    # int *cMask (full ny*nx mask array)
    ]
    convolve_c.restype = None

    ny, nx = image.shape

    # Convert to C pointers
    img_ptr = array_to_cptr(image)
    out_ptr = array_to_cptr(output)
    coeffs_ptr = array_to_double_cptr(kernel_coeffs)

    # spatial_convolve takes float** variance and immediately dereferences it
    # (if (*variance) == NULL → no variance image). We must pass a pointer to
    # a NULL float* — not NULL itself — to avoid a segfault on *variance.
    null_var_ptr = ctypes.c_void_p(None)

    # Set up globals needed by spatial_convolve: rPixX, rPixY, mRData, kcStep.
    # setupSpatialConvolve allocates a zeroed mask buffer, sets mRData to it,
    # and sets rPixX/rPixY so make_kernel_local can normalise pixel coordinates.
    # C convention: xSize=nx (width), ySize=ny (height); setupSpatialConvolve
    # uses the same convention (rPixX=nx, rPixY=ny).
    mask_ptr = lib.setupSpatialConvolve(nx, ny)
    try:
        convolve_c(img_ptr, ctypes.byref(null_var_ptr), nx, ny, coeffs_ptr, out_ptr, mask_ptr)
    finally:
        lib.cleanupSpatialConvolve()

    return output
