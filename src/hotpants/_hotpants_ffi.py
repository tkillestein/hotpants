"""
HOTPANTS cffi FFI Declarations

This module defines the C API interface using cffi's out-of-line mode.
It parses hotpants_api.h and generates compiled extension modules.

Note: This file is processed by cffi during package build. The cffi library
compiles it into a native extension (_hotpants_cffi.so) that provides fast
access to the C code without overhead.
"""

from cffi import FFI

ffi = FFI()

# Define C types and functions from the public API header
ffi.cdef("""
    /* Stamp data structure */
    typedef struct {
        int       x0, y0;
        int       x, y;
        int       nx, ny;
        int       *xss;
        int       *yss;
        int       nss;
        int       sscnt;
        double    **vectors;
        double    *krefArea;
        double    **mat;
        double    *scprod;
        double    sum;
        double    mean;
        double    median;
        double    mode;
        double    sd;
        double    fwhm;
        double    lfwhm;
        double    chi2;
        double    norm;
        double    diff;
    } stamp_struct;

    /* Core API functions */
    int allocateStamps(stamp_struct *stamps, int n);
    void freeStampMem(stamp_struct *stamps, int n);

    void buildStamps(int sXMin, int sXMax, int sYMin, int sYMax,
                     int *rPixX, int *rPixY, int nStampX, int nStampY,
                     int hwKSStamp,
                     stamp_struct *stamps, stamp_struct *tStamps, stamp_struct *iStamps,
                     float *iData, float *tData,
                     float tUThresh, float tLThresh);

    int getPsfCenters(stamp_struct *stamp, float *iData, int nsx, int nsy,
                      double smin, int nKSStamps, int nKSEdge);

    int getStampStats3(float *data, int nx, int ny, int nsy,
                       int stat_type,
                       double *mean, double *median, double *mode,
                       double *sd, double *fwhm, double *lfwhm,
                       int verbose, int nThread, int nComp);

    void fitKernel(stamp_struct *stamps, float *imRef, float *imConv,
                   float *imNoise, double *kernel_coeffs,
                   double *meansig, double *scatter, int *n_skipped);

    void spatial_convolve(float *image, float **var_image, int ny, int nx,
                          double *kernel_coeffs, float *output, int *conv_method);

    /* Global configuration variables */
    extern int       hwKernel;
    extern int       kerOrder;
    extern int       bgOrder;
    extern float     kerFitThresh;
    extern float     scaleFitThresh;
    extern int       nKSStamps;
    extern int       hwKSStamp;
    extern int       nStampX, nStampY;
    extern int       useFullSS;

    extern float     tUThresh, tLThresh;
    extern float     iUThresh, iLThresh;
    extern float     tGain, iGain;
    extern float     tRdnoise, iRdnoise;
    extern float     tPedestal, iPedestal;

    extern int       nCompKer;
    extern int       nComp;
    extern int       nCompBG;

    extern int       verbose;
    extern int       nThread;

    /* Memory allocation helpers (for creating double pointers, etc.) */
    double** malloc_double_array(int rows, int cols);
    void free_double_array(double** arr, int rows);
""")

# Configure source code to compile
# Set source to compile the C code with cffi
ffi.set_source(
    "hotpants._hotpants_cffi",  # module name
    """
    #define HOTPANTS_DEFINE_GLOBALS
    #include "hotpants_api.h"

    /* Simple allocators for double pointers (needed for Python to pass to C) */
    double** malloc_double_array(int rows, int cols) {
        double** arr = (double**)malloc(rows * sizeof(double*));
        for (int i = 0; i < rows; i++) {
            arr[i] = (double*)malloc(cols * sizeof(double));
        }
        return arr;
    }

    void free_double_array(double** arr, int rows) {
        for (int i = 0; i < rows; i++) {
            free(arr[i]);
        }
        free(arr);
    }
    """,
    sources=[
        "src/alard.c",
        "src/functions.c",
    ],
    include_dirs=["src"],
    libraries=["cfitsio", "lapack", "lapacke", "m"],
    extra_compile_args=["-O3", "-march=native", "-std=c17"],
    language="c",
)

if __name__ == "__main__":
    ffi.compile(verbose=True)
