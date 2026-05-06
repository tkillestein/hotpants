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
import os
from functools import lru_cache
from pathlib import Path


# TODO: resolve C901 error
def load_hotpants_library() -> ctypes.CDLL:  # noqa: C901
    """
    Load the pre-compiled HOTPANTS C library.

    Searches in standard locations and environment paths.
    Returns the ctypes CDLL instance.

    Raises:
        OSError: if library not found
    """

    lib_name = "libhotpants"

    # Try to find libhotpants using ctypes' find_library (searches LD_LIBRARY_PATH, etc.)
    lib_path = ctypes.util.find_library(lib_name)
    if lib_path:
        try:
            return ctypes.CDLL(lib_path)
        except OSError:
            pass  # Fall through to manual search

    # Try build directory relative to this file
    possible_paths = [
        Path(__file__).parent.parent.parent / "build",
        Path(__file__).parent.parent.parent / "build" / "src",
        Path("/usr/local/lib"),
        Path("/usr/lib"),
    ]

    # Also check LD_LIBRARY_PATH
    ld_path = os.environ.get("LD_LIBRARY_PATH", "")
    if ld_path:
        for path_str in ld_path.split(":"):
            if path_str:
                possible_paths.insert(0, Path(path_str))

    for path_dir in possible_paths:
        if path_dir.exists():
            for pattern in [f"{lib_name}.so*", f"{lib_name}.dylib", f"{lib_name}.dll"]:
                matches = list(path_dir.glob(pattern))
                if matches:
                    lib_file = matches[0]
                    try:
                        return ctypes.CDLL(str(lib_file))
                    except OSError:
                        # Continue to next path if this one doesn't load
                        continue

    # Give helpful error message
    searched_paths = "\n  ".join(str(p) for p in possible_paths)
    msg = (
        f"Could not find {lib_name}. Searched:\n  {searched_paths}\n"
        "Please build the C library first:\n"
        "  cmake -B build && cmake --build build\n"
        "Then set LD_LIBRARY_PATH (or run tests with):\n"
        "  export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH"
    )
    raise OSError(msg)


@lru_cache(maxsize=1)
def get_library() -> ctypes.CDLL:
    return load_hotpants_library()


# =====================================================================
# ctypes Wrapper Functions
# =====================================================================


def get_hotpants_library_functions() -> dict[str, object]:
    """
    Define C function signatures for ctypes.

    Returns a dict mapping function names to ctypes function wrappers
    with proper argument and return types.
    """
    lib = get_library()

    functions = {}

    # int allocateStamps(stamp_struct *stamps, int n)
    functions["allocateStamps"] = lib.allocateStamps
    functions["allocateStamps"].argtypes = [ctypes.c_void_p, ctypes.c_int]
    functions["allocateStamps"].restype = ctypes.c_int

    # void freeStampMem(stamp_struct *stamps, int n)
    functions["freeStampMem"] = lib.freeStampMem
    functions["freeStampMem"].argtypes = [ctypes.c_void_p, ctypes.c_int]
    functions["freeStampMem"].restype = None

    # void buildStamps(int sXMin, sXMax, sYMin, sYMax, int *rPixX, *rPixY,
    #                  int nStampX, nStampY, hwKSStamp,
    #                  stamp_struct *stamps, *tStamps, *iStamps,
    #                  float *iData, *tData, tUThresh, tLThresh)
    functions["buildStamps"] = lib.buildStamps
    functions["buildStamps"].argtypes = [
        ctypes.c_int,  # sXMin
        ctypes.c_int,  # sXMax
        ctypes.c_int,  # sYMin
        ctypes.c_int,  # sYMax
        ctypes.c_void_p,  # int *rPixX
        ctypes.c_void_p,  # int *rPixY
        ctypes.c_int,  # nStampX
        ctypes.c_int,  # nStampY
        ctypes.c_int,  # hwKSStamp
        ctypes.c_void_p,  # stamp_struct *stamps
        ctypes.c_void_p,  # stamp_struct *tStamps
        ctypes.c_void_p,  # stamp_struct *iStamps
        ctypes.c_void_p,  # float *iData
        ctypes.c_void_p,  # float *tData
        ctypes.c_float,  # tUThresh
        ctypes.c_float,  # tLThresh
    ]
    functions["buildStamps"].restype = None

    # int getPsfCenters(stamp_struct *stamp, float *iData, int nsx, nsy,
    #                   double smin, int nKSStamps, nKSEdge)
    functions["getPsfCenters"] = lib.getPsfCenters
    functions["getPsfCenters"].argtypes = [
        ctypes.c_void_p,  # stamp_struct *stamp
        ctypes.c_void_p,  # float *iData
        ctypes.c_int,  # nsx
        ctypes.c_int,  # nsy
        ctypes.c_double,  # smin
        ctypes.c_int,  # nKSStamps
        ctypes.c_int,  # nKSEdge
    ]
    functions["getPsfCenters"].restype = ctypes.c_int

    # int getStampStats3(...)
    functions["getStampStats3"] = lib.getStampStats3
    functions["getStampStats3"].argtypes = [
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
    functions["getStampStats3"].restype = ctypes.c_int

    # void fitKernel(stamp_struct *stamps, float *imRef, *imConv, *imNoise,
    #                double *kernel_coeffs, *meansig, *scatter, int *n_skipped)
    functions["fitKernel"] = lib.fitKernel
    functions["fitKernel"].argtypes = [
        ctypes.c_void_p,  # stamp_struct *stamps
        ctypes.c_void_p,  # float *imRef
        ctypes.c_void_p,  # float *imConv
        ctypes.c_void_p,  # float *imNoise (nullable)
        ctypes.c_void_p,  # double *kernel_coeffs
        ctypes.POINTER(ctypes.c_double),  # *meansig
        ctypes.POINTER(ctypes.c_double),  # *scatter
        ctypes.POINTER(ctypes.c_int),  # *n_skipped
    ]
    functions["fitKernel"].restype = None

    # void spatial_convolve(float *image, float **var_image, int ny, nx,
    #                       double *kernel_coeffs, float *output, int *conv_method)
    functions["spatial_convolve"] = lib.spatial_convolve
    functions["spatial_convolve"].argtypes = [
        ctypes.c_void_p,  # float *image
        ctypes.c_void_p,  # float **var_image (nullable)
        ctypes.c_int,  # ny
        ctypes.c_int,  # nx
        ctypes.c_void_p,  # double *kernel_coeffs
        ctypes.c_void_p,  # float *output
        ctypes.POINTER(ctypes.c_int),  # *conv_method
    ]
    functions["spatial_convolve"].restype = None

    # Wrapper functions for buildStamps context management
    # int initBuildStampsContext(float* template, float* science,
    #                            int ny, int nx, int n_regions_x, int n_regions_y,
    #                            int stamps_per_region_x, int stamps_per_region_y)
    functions["initBuildStampsContext"] = lib.initBuildStampsContext
    functions["initBuildStampsContext"].argtypes = [
        ctypes.c_void_p,  # float *template
        ctypes.c_void_p,  # float *science
        ctypes.c_int,  # ny
        ctypes.c_int,  # nx
        ctypes.c_int,  # n_regions_x
        ctypes.c_int,  # n_regions_y
        ctypes.c_int,  # stamps_per_region_x
        ctypes.c_int,  # stamps_per_region_y
    ]
    functions["initBuildStampsContext"].restype = ctypes.c_int

    # int buildStampsRegion(int region_x, int region_y,
    #                       stamp_struct* stamps, int n_stamps)
    functions["buildStampsRegion"] = lib.buildStampsRegion
    functions["buildStampsRegion"].argtypes = [
        ctypes.c_int,  # region_x
        ctypes.c_int,  # region_y
        ctypes.c_void_p,  # stamp_struct *stamps
        ctypes.c_int,  # n_stamps
    ]
    functions["buildStampsRegion"].restype = ctypes.c_int

    # void cleanupBuildStampsContext(void)
    functions["cleanupBuildStampsContext"] = lib.cleanupBuildStampsContext
    functions["cleanupBuildStampsContext"].argtypes = []
    functions["cleanupBuildStampsContext"].restype = None

    return functions


# Global variable access via ctypes
def get_global_int(name: str) -> int:
    """Get value of integer global variable."""
    lib = get_library()
    var = ctypes.c_int.in_dll(lib, name)
    return var.value


def set_global_int(name: str, value: int) -> None:
    """Set value of integer global variable."""
    lib = get_library()
    var = ctypes.c_int.in_dll(lib, name)
    var.value = value


def get_global_float(name: str) -> float:
    """Get value of float global variable."""
    lib = get_library()
    var = ctypes.c_float.in_dll(lib, name)
    return var.value


def set_global_float(name: str, value: float) -> None:
    """Set value of float global variable."""
    lib = get_library()
    var = ctypes.c_float.in_dll(lib, name)
    var.value = value
