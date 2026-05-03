"""
HOTPANTS ctypes FFI Layer

This module provides a ctypes interface to the HOTPANTS C library.
No compilation is needed — ctypes loads the pre-compiled C library at runtime.

The C library must be built separately via CMake:
    cmake -B build && cmake --build build

Then loaded by finding the shared library (.so / .dll / .dylib).
"""

import ctypes
import ctypes.util
import sys
from pathlib import Path


def load_hotpants_library():
    """
    Load the pre-compiled HOTPANTS C library.

    Searches in standard locations and environment paths.
    Returns the ctypes CDLL instance.

    Raises:
        OSError: if library not found
    """
    # Try to find libhotpants (built by CMake)
    lib_name = "libhotpants"

    # Try standard locations first
    lib_path = ctypes.util.find_library(lib_name)
    if lib_path:
        return ctypes.CDLL(lib_path)

    # Try build directory
    possible_paths = [
        Path(__file__).parent.parent.parent / "build",
        Path(__file__).parent.parent.parent / "build" / "src",
        Path("/usr/local/lib"),
        Path("/usr/lib"),
    ]

    for path_dir in possible_paths:
        if path_dir.exists():
            for pattern in [f"{lib_name}.so*", f"{lib_name}.dylib", f"{lib_name}.dll"]:
                matches = list(path_dir.glob(pattern))
                if matches:
                    return ctypes.CDLL(str(matches[0]))

    raise OSError(
        f"Could not find {lib_name}. "
        "Please build the C library first:\n"
        "  cmake -B build && cmake --build build\n"
        "Then set LD_LIBRARY_PATH to the build directory."
    )


# Load the library (deferred until first use)
_lib = None


def get_library():
    """Get the loaded HOTPANTS library, loading it if necessary."""
    global _lib
    if _lib is None:
        _lib = load_hotpants_library()
    return _lib


# =====================================================================
# ctypes Wrapper Functions
# =====================================================================

def get_hotpants_library_functions():
    """
    Define C function signatures for ctypes.

    Returns a dict mapping function names to ctypes function wrappers
    with proper argument and return types.
    """
    lib = get_library()

    functions = {}

    # int allocateStamps(stamp_struct *stamps, int n)
    functions['allocateStamps'] = lib.allocateStamps
    functions['allocateStamps'].argtypes = [ctypes.c_void_p, ctypes.c_int]
    functions['allocateStamps'].restype = ctypes.c_int

    # void freeStampMem(stamp_struct *stamps, int n)
    functions['freeStampMem'] = lib.freeStampMem
    functions['freeStampMem'].argtypes = [ctypes.c_void_p, ctypes.c_int]
    functions['freeStampMem'].restype = None

    # int getStampStats3(...)
    functions['getStampStats3'] = lib.getStampStats3
    functions['getStampStats3'].argtypes = [
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
    functions['getStampStats3'].restype = ctypes.c_int

    return functions


# Global variable access via ctypes
def get_global_int(name):
    """Get value of integer global variable."""
    lib = get_library()
    var = ctypes.c_int.in_dll(lib, name)
    return var.value


def set_global_int(name, value):
    """Set value of integer global variable."""
    lib = get_library()
    var = ctypes.c_int.in_dll(lib, name)
    var.value = value


def get_global_float(name):
    """Get value of float global variable."""
    lib = get_library()
    var = ctypes.c_float.in_dll(lib, name)
    return var.value


def set_global_float(name, value):
    """Set value of float global variable."""
    lib = get_library()
    var = ctypes.c_float.in_dll(lib, name)
    var.value = value
