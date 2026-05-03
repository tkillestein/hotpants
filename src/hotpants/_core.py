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

import numpy as np
import ctypes
from contextlib import contextmanager
from . import _hotpants_ffi


def validate_image_array(arr, name="image"):
    """
    Validate image array for use with C code.

    Checks:
    - dtype is float32
    - array is 2D
    - array is C-contiguous
    - no NaN or Inf values

    Args:
        arr: numpy array to validate
        name: name for error messages

    Raises:
        TypeError: if dtype is not float32
        ValueError: if array shape or values are invalid
    """
    if arr.dtype != np.float32:
        raise TypeError(
            f"{name}: expected float32, got {arr.dtype}. "
            f"Convert with arr.astype(np.float32)."
        )
    if arr.ndim != 2:
        raise ValueError(
            f"{name}: expected 2D array, got shape {arr.shape}"
        )
    if not arr.flags['C_CONTIGUOUS']:
        raise ValueError(
            f"{name}: array must be C-contiguous (row-major). "
            "Use np.ascontiguousarray(arr) to convert."
        )
    if not np.isfinite(arr).all():
        raise ValueError(
            f"{name}: contains NaN or Inf values. "
            "Check for bad pixels or invalid thresholds."
        )


def array_to_cptr(arr):
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
        raise TypeError(f"Expected float32, got {arr.dtype}")
    if not arr.flags['C_CONTIGUOUS']:
        arr = np.ascontiguousarray(arr)

    # Cast numpy array pointer to ctypes void*
    # This is a borrowed reference; numpy array owns the memory
    return ctypes.cast(arr.ctypes.data, ctypes.c_void_p)


def array_to_double_cptr(arr):
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
        raise TypeError(f"Expected float32/float64, got {arr.dtype}")

    if not arr.flags['C_CONTIGUOUS']:
        arr = np.ascontiguousarray(arr)

    return ctypes.cast(arr.ctypes.data, ctypes.c_void_p)


def allocate_array(shape, dtype=np.float32):
    """
    Allocate a numpy array suitable for C code.

    Arrays are always C-contiguous (row-major), required by C code.

    Args:
        shape: tuple (ny, nx)
        dtype: data type (default: float32)

    Returns:
        numpy array, C-contiguous, zero-initialized
    """
    return np.zeros(shape, dtype=dtype, order='C')


@contextmanager
def global_state(config_dict):
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
        'hwKernel', 'kerOrder', 'bgOrder', 'nKSStamps', 'hwKSStamp',
        'nRegX', 'nRegY', 'nStampX', 'nStampY', 'useFullSS',
        'verbose', 'nThread', 'nCompKer', 'nComp', 'nCompBG',
    }
    float_vars = {
        'kerFitThresh', 'scaleFitThresh',
        'tUThresh', 'tLThresh', 'iUThresh', 'iLThresh',
        'tGain', 'iGain', 'tRdnoise', 'iRdnoise', 'tPedestal', 'iPedestal',
    }

    # Save original values
    saved = {}
    for key, val in config_dict.items():
        if key in int_vars:
            saved[key] = _hotpants_ffi.get_global_int(key)
        elif key in float_vars:
            saved[key] = _hotpants_ffi.get_global_float(key)
        else:
            raise AttributeError(f"Unknown global variable: {key}")

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
# Stamp Management Wrappers
# =====================================================================

def allocate_stamps(n_stamps):
    """
    Allocate stamp array via ctypes.

    Args:
        n_stamps: number of stamps

    Returns:
        ctypes array of stamp structures (opaque to Python)

    Raises:
        RuntimeError: if allocation fails
    """
    # For now, return a placeholder
    # Full implementation requires stamp_struct definition in ctypes
    # This will be implemented when C integration is complete
    raise NotImplementedError("Stamp allocation requires C library binding")


def free_stamps(stamps, n_stamps):
    """Free stamp array memory."""
    # Placeholder for C integration
    pass


# =====================================================================
# Image Statistics Wrapper
# =====================================================================

def get_stamp_stats(data_array, verbose=0, n_thread=1, n_comp=3):
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
    getStampStats3 = lib.getStampStats3
    getStampStats3.argtypes = [
        ctypes.c_void_p,  # float *data
        ctypes.c_int,     # nx
        ctypes.c_int,     # ny
        ctypes.c_int,     # nsy
        ctypes.c_int,     # stat_type
        ctypes.POINTER(ctypes.c_double),  # *mean
        ctypes.POINTER(ctypes.c_double),  # *median
        ctypes.POINTER(ctypes.c_double),  # *mode
        ctypes.POINTER(ctypes.c_double),  # *sd
        ctypes.POINTER(ctypes.c_double),  # *fwhm
        ctypes.POINTER(ctypes.c_double),  # *lfwhm
        ctypes.c_int,     # verbose
        ctypes.c_int,     # nThread
        ctypes.c_int,     # nComp
    ]
    getStampStats3.restype = ctypes.c_int

    # Call C function
    n_used = getStampStats3(
        data_ptr, nx, ny, ny,  # nx, ny, nsy (row count)
        0,  # stat_type (unused)
        ctypes.byref(mean_val),
        ctypes.byref(median_val),
        ctypes.byref(mode_val),
        ctypes.byref(sd_val),
        ctypes.byref(fwhm_val),
        ctypes.byref(lfwhm_val),
        verbose, n_thread, n_comp
    )

    return {
        'n_used': n_used,
        'mean': mean_val.value,
        'median': median_val.value,
        'mode': mode_val.value,
        'sd': sd_val.value,
        'fwhm': fwhm_val.value,
        'lfwhm': lfwhm_val.value,
    }


# =====================================================================
# Polynomial Basis Term Count
# =====================================================================

def num_polynomial_terms(order):
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


def get_kernel_info():
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
        'kernel_components': _hotpants_ffi.get_global_int('nCompKer'),
        'kernel_order': _hotpants_ffi.get_global_int('kerOrder'),
        'bg_order': _hotpants_ffi.get_global_int('bgOrder'),
        'kernel_terms': _hotpants_ffi.get_global_int('nComp'),
        'bg_terms': _hotpants_ffi.get_global_int('nCompBG'),
    }
