"""
Low-level cffi wrappers and numpy-to-C conversions.

This module provides:
- Array conversion utilities (numpy ↔ C pointers)
- Global state management (set/restore configuration variables)
- Low-level wrapper functions around core C functions
- Error handling and validation

Functions here are not intended for direct use by Python users;
they are internal utilities called by the high-level API in api.py.
"""

import numpy as np
from contextlib import contextmanager

# Import the compiled cffi module
try:
    from hotpants import _hotpants_cffi as lib
    ffi = lib.ffi
except ImportError as e:
    raise ImportError(
        "HOTPANTS cffi extension not found. "
        "Did you forget to install? Try: pip install -e ."
    ) from e


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
    Convert numpy array to cffi float pointer (borrowed reference).

    Creates a C float* that points to the numpy array's data.
    The numpy array must stay in scope during the C function call.

    Args:
        arr: numpy array (float32, 2D, C-contiguous)

    Returns:
        cffi float* pointer to array data
    """
    if arr.dtype != np.float32:
        raise TypeError(f"Expected float32, got {arr.dtype}")
    if not arr.flags['C_CONTIGUOUS']:
        arr = np.ascontiguousarray(arr)

    # Cast numpy array pointer to cffi float*
    # This is a borrowed reference; numpy array owns the memory
    return ffi.cast("float *", arr.ctypes.data)


def array_to_double_cptr(arr):
    """
    Convert numpy array (float32 or float64) to cffi double pointer.

    If input is float32, creates a temporary float64 copy.
    The resulting pointer is valid only during C function execution.

    Args:
        arr: numpy array (float32 or float64, 2D, C-contiguous)

    Returns:
        cffi double* pointer
    """
    if arr.dtype == np.float32:
        arr = arr.astype(np.float64)
    elif arr.dtype != np.float64:
        raise TypeError(f"Expected float32/float64, got {arr.dtype}")

    if not arr.flags['C_CONTIGUOUS']:
        arr = np.ascontiguousarray(arr)

    return ffi.cast("double *", arr.ctypes.data)


def cptr_to_array(ptr, shape):
    """
    Wrap C pointer as numpy array (borrowed memory, read-only view).

    Creates a numpy array view over C-allocated memory. The view is
    read-only to prevent accidental modification. The C code must
    manage the lifetime of the underlying buffer.

    Args:
        ptr: cffi float* or double* pointer
        shape: tuple (ny, nx)

    Returns:
        numpy array with shape (ny, nx), read-only view
    """
    ny, nx = shape
    dtype = np.float32 if 'float' in str(ffi.typeof(ptr)) else np.float64
    size = ny * nx

    # Create buffer view from C pointer
    buffer_view = ffi.buffer(ptr, size * (4 if dtype == np.float32 else 8))
    arr = np.frombuffer(buffer_view, dtype=dtype)

    # Reshape and return read-only view
    arr = arr.reshape(shape)
    arr.flags.writeable = False
    return arr


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
    # Save original values
    saved = {}
    for key in config_dict:
        if not hasattr(lib, key):
            raise AttributeError(f"Unknown global variable: {key}")
        saved[key] = getattr(lib, key)

    # Set new values
    for key, val in config_dict.items():
        setattr(lib, key, val)

    try:
        yield
    finally:
        # Restore original values
        for key, val in saved.items():
            setattr(lib, key, val)


# =====================================================================
# Stamp Management Wrappers
# =====================================================================

def allocate_stamps(n_stamps):
    """
    Allocate stamp array via cffi.

    Args:
        n_stamps: number of stamps

    Returns:
        cffi stamp_struct* array

    Raises:
        RuntimeError: if allocation fails
    """
    stamps = ffi.new("stamp_struct[]", n_stamps)
    result = lib.allocateStamps(stamps, n_stamps)
    if result != 0:
        raise RuntimeError(f"allocateStamps failed with code {result}")
    return stamps


def free_stamps(stamps, n_stamps):
    """Free stamp array memory."""
    lib.freeStampMem(stamps, n_stamps)


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

    # Allocate output variables
    mean_ptr = ffi.new("double *")
    median_ptr = ffi.new("double *")
    mode_ptr = ffi.new("double *")
    sd_ptr = ffi.new("double *")
    fwhm_ptr = ffi.new("double *")
    lfwhm_ptr = ffi.new("double *")

    # Call C function
    n_used = lib.getStampStats3(
        data_ptr, nx, ny, ny,  # nx, ny, nsy (row count)
        0,  # stat_type (unused)
        mean_ptr, median_ptr, mode_ptr, sd_ptr, fwhm_ptr, lfwhm_ptr,
        verbose, n_thread, n_comp
    )

    return {
        'n_used': n_used,
        'mean': mean_ptr[0],
        'median': median_ptr[0],
        'mode': mode_ptr[0],
        'sd': sd_ptr[0],
        'fwhm': fwhm_ptr[0],
        'lfwhm': lfwhm_ptr[0],
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
        'kernel_components': lib.nCompKer,
        'kernel_order': lib.kerOrder,
        'bg_order': lib.bgOrder,
        'kernel_terms': lib.nComp,
        'bg_terms': lib.nCompBG,
    }
