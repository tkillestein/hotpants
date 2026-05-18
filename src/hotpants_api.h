#ifndef HOTPANTS_API_H
#define HOTPANTS_API_H

/*
 * HOTPANTS Python API — Minimal C Interface
 *
 * This header exposes the core kernel-fitting and convolution functions
 * for Python binding via cffi. It includes only public data structures
 * and function signatures necessary for the Python API.
 *
 * For implementation details and additional functions, see:
 * - alard.c: kernel fitting algorithm (Alard & Lupton 1998)
 * - functions.c: stamp utilities, statistics, convolution
 * - globals.h: global variable declarations
 */

#include <stddef.h>

/* =====================================================================
 * Core Data Structure: Stamp
 * =====================================================================
 * A "stamp" represents a small cutout of the image (typically ~32x32 px)
 * within a region. The kernel is fit separately in each stamp, and the
 * results are combined to form the spatially-varying kernel.
 *
 * Reference: Alard & Lupton (1998), Section 2.1 (kernel parameterization).
 */
typedef struct {
  int x0, y0;       /* origin of stamp in region coords */
  int x, y;         /* center of stamp in region coords */
  int nx, ny;       /* size of stamp in pixels */
  int* xss;         /* x location of kernel test substamp centers */
  int* yss;         /* y location of kernel test substamp centers */
  int nss;          /* number of detected kernel substamps (1 .. nss) */
  int sscnt;        /* which substamp to use (0 .. nss-1) */
  double** vectors; /* convolved template basis vectors for this stamp */
  double* krefArea; /* reference (template) kernel substamp data */
  double** mat;     /* fitting matrix (normal equations) */
  double* scprod;   /* scalar products (right-hand side) */
  double sum;       /* sum of fabs(data), used for sigma scaling */
  double mean;      /* stamp mean value */
  double median;    /* stamp median value */
  double mode;      /* estimated sky mode (background) */
  double sd;        /* standard deviation */
  double fwhm;      /* full width half max of PSF */
  double lfwhm;     /* log(fwhm) */
  double chi2;      /* chi-squared residual from kernel fit */
  double norm;      /* kernel sum (integral) */
  double diff;      /* (norm - mean_ksum) * sqrt(sum) */
} stamp_struct;

/* =====================================================================
 * Configuration Struct — Thread-Safe API
 * =====================================================================
 * Bundles all configuration inputs for a single fit_kernel / spatial_convolve
 * call. Pass a pointer to the _ctx wrapper functions below. Computed fields
 * (nCompKer, nComp, ...) are written back by initKernelGlobals_ctx().
 *
 * Using this struct avoids setting individual C globals from Python, enabling
 * safe concurrent calls via the internal mutex.
 */
typedef struct {
    /* Kernel geometry */
    int    hwKernel;
    int    kerOrder;
    int    bgOrder;

    /* Basis type */
    int    iBasisType;
    double rDeltaRegularization;
    int    useTPS;
    double tpsSmoothing;

    /* Stamp layout */
    int    nKSStamps;
    int    hwKSStamp;
    int    nRegX;
    int    nRegY;
    int    nStampX;
    int    nStampY;
    int    useFullSS;

    /* Image thresholds */
    float  tUThresh;
    float  tLThresh;
    float  iUThresh;
    float  iLThresh;

    /* Noise model */
    float  tGain;
    float  iGain;
    float  tRdnoise;
    float  iRdnoise;
    float  tPedestal;
    float  iPedestal;

    /* Fit quality */
    float  kerFitThresh;
    float  scaleFitThresh;

    /* Execution control */
    int    verbose;
    int    nThread;

    /* Derived fields written by initKernelGlobals_ctx() */
    int    nCompKer;
    int    nComp;
    int    nCompBG;
    int    nCompTotal;
    int    fwKernel;
    int    fwKSStamp;
    int    fwStamp;
    int    kcStep;

    /* Valid stamp count: set by build phase, read by fit phase */
    int    nS;
} hotpants_config_t;

/* =====================================================================
 * Core API Functions
 * =====================================================================
 * These functions comprise the minimal public interface for kernel fitting
 * and convolution. All other functions in alard.c/functions.c are treated
 * as internal implementation details.
 */

/*
 * Allocate and initialize stamp array.
 *
 * Args:
 *   stamps: pointer to stamp array of size n
 *   n: number of stamps to allocate
 *
 * Returns:
 *   0 on success, -1 on allocation failure
 */
int allocateStamps(stamp_struct* stamps, int n);

/*
 * Free internal memory associated with stamp array.
 *
 * Args:
 *   stamps: stamp array to free
 *   n: number of stamps
 */
void freeStampMem(stamp_struct* stamps, int n);

/*
 * Divide image region into stamp grid and identify bright PSF centers.
 *
 * This function:
 * 1. Lays out a grid of stamps (nStampX x nStampY per region)
 * 2. Within each stamp, identifies nKSStamps bright star centers
 * 3. Populates stamp_struct with pixel locations and statistics
 *
 * Args:
 *   sXMin, sXMax: X pixel range of region in full image
 *   sYMin, sYMax: Y pixel range of region in full image
 *   iData: science image array (ny, nx) in pixels
 *   tData: template image array (ny, nx) in pixels
 *   tUThresh, tLThresh: template upper/lower valid data thresholds
 *   tGain: template gain (e-/ADU)
 *   tRdnoise: template readnoise (e-)
 *   stamps: output array of stamp_struct, already allocated
 *   tStamps, iStamps: auxiliary stamp arrays
 *
 * Reference: Alard & Lupton (1998), Section 2.1
 */
void buildStamps(int sXMin, int sXMax, int sYMin, int sYMax, int* rPixX,
                 int* rPixY, int nStampX, int nStampY, int hwKSStamp,
                 stamp_struct* stamps, stamp_struct* tStamps,
                 stamp_struct* iStamps, float* iData, float* tData,
                 float tUThresh, float tLThresh);

/*
 * Locate bright PSF centers (substamp locations) within a stamp.
 *
 * This function identifies nKSStamps bright point sources (stars) within
 * a single stamp. These locations are used to constrain the kernel fit.
 *
 * Args:
 *   stamp: single stamp_struct to process
 *   iData: science image array
 *   nsx, nsy: image dimensions
 *   smin: minimum significant value for source detection
 *   nKSStamps: number of substamps to locate
 *   nKSEdge: edge buffer to avoid border effects
 *
 * Returns:
 *   Number of substamps successfully located (may be < nKSStamps)
 *
 * Reference: Alard & Lupton (1998), Section 2.1
 */
int getPsfCenters(stamp_struct* stamp, float* iData, int nsx, int nsy,
                  double smin, int nKSStamps, int nKSEdge);

/*
 * Compute robust image statistics (mean, median, mode, FWHM).
 *
 * This function estimates image statistics using histogram analysis with
 * sigma-clipping to handle outliers robustly. Returns FWHM estimate
 * suitable for PSF characterization.
 *
 * Args:
 *   data: pixel array
 *   x0Reg, y0Reg: region origin (for internal reference)
 *   nPixX, nPixY: array dimensions (width, height)
 *   sum, mean, median, mode, sd, fwhm, lfwhm: output statistics
 *   umask: pixel mask upper threshold (e.g., 0xff for OK pixels)
 *   smask: pixel mask significance threshold (e.g., 0x8000)
 *   maxiter: maximum iterations for sigma-clipping
 *
 * Returns:
 *   Number of pixels used in final estimate (after sigma-clipping)
 *
 * Reference: functions.c getStampStats3() for histogram algorithm.
 */
int getStampStats3(float* data, int x0Reg, int y0Reg, int nPixX, int nPixY,
                   double* sum, double* mean, double* median, double* mode,
                   double* sd, double* fwhm, double* lfwhm, int umask,
                   int smask, int maxiter);

/*
 * Fit spatially-varying convolution kernel via least-squares.
 *
 * This function solves the normal equations derived from Alard & Lupton's
 * algorithm. It processes all stamps in a region, accumulates the normal
 * equation matrix, performs Cholesky decomposition (LAPACK), and returns
 * the fitted kernel coefficients.
 *
 * The kernel is parameterized as:
 *   K(x,y) = Σᵢ cᵢ(x,y) * φᵢ(x,y)
 * where φᵢ are fixed Gaussian basis functions and cᵢ are spatial polynomials.
 *
 * Args:
 *   stamps: array of stamp_struct with convolution vectors pre-computed
 *   imRef: reference (template) image array
 *   imConv: image to be convolved (usually science image)
 *   imNoise: optional noise image (NULL if not used)
 *   kernel_coeffs: output array of fitted kernel coefficients
 *   meansig: output mean significance of fit
 *   scatter: output scatter in fit significance
 *   n_skipped: output number of stamps excluded due to poor fit
 *
 * Returns:
 *   (void; check n_skipped and return values via pointers)
 *
 * Reference: Alard & Lupton (1998), Section 2.2 (kernel fitting)
 */
void fitKernel(stamp_struct* stamps, float* imRef, float* imConv,
               float* imNoise, double* kernel_coeffs, double* meansig,
               double* scatter, int* n_skipped);

/*
 * Apply spatially-varying kernel as convolution to full image.
 *
 * This function evaluates the kernel at every pixel position and performs
 * the convolution I ⊗ K(x,y), storing the result in output array.
 *
 * The kernel is reconstructed from coefficients:
 *   K(x,y) = Σᵢ cᵢ(x,y) * φᵢ
 * Each basis element φᵢ is convolved with the input image via FFT or
 * direct convolution (configurable at runtime).
 *
 * Args:
 *   image: input image array (ny, nx)
 *   var_image: per-pixel variance array (or NULL)
 *   ny, nx: image dimensions
 *   kernel_coeffs: fitted kernel coefficients from fitKernel()
 *   output: pre-allocated output array (ny, nx), filled in-place
 *   conv_method: convolution method flag (0=direct, 1=FFT, etc.)
 *
 * Returns:
 *   (void)
 *
 * Reference: Alard & Lupton (1998), Section 3 (convolution)
 */
void spatial_convolve(float* image, float** var_image, int ny, int nx,
                      double* kernel_coeffs, float* output, int* conv_method);

/* =====================================================================
 * Wrapper Functions for Python Integration
 * =====================================================================
 * Two families of wrappers are provided:
 *
 * 1. Legacy (original) functions — still available; set globals directly.
 *    Used by main.c CLI and any code that sets globals before calling.
 *
 * 2. _ctx variants — accept a hotpants_config_t* as the first argument.
 *    They apply config → globals atomically under an internal mutex, call
 *    the algorithm, and read back computed outputs into the struct.
 *    Prefer these for Python callers; safe for concurrent threads.
 */

/*
 * Initialize derived kernel globals from hwKernel/kerOrder/bgOrder.
 *
 * Must be called before allocateStamps() or buildStamps(). Computes:
 *   nCompKer, nComp, nC, fwKSStamp, fwKernel, fwStamp
 * and initialises the Gaussian basis arrays (ngauss, deg_fixe, sigma_gauss)
 * with standard defaults matching the 3-Gaussian HOTPANTS basis.
 *
 * The hwKernel, kerOrder, bgOrder, nKSStamps, hwKSStamp globals must be set
 * by the caller before invoking this function.
 *
 * Args:
 *   image_nx, image_ny: full image dimensions
 *   n_reg_x, n_reg_y: number of regions per axis
 *   n_stamp_x, n_stamp_y: stamps per region per axis
 *
 * Returns:
 *   0 on success, -1 on allocation failure
 */
int initKernelGlobals(int image_nx, int image_ny,
                      int n_reg_x, int n_reg_y,
                      int n_stamp_x, int n_stamp_y);

/*
 * Initialize buildStamps context for a complete image.
 *
 * Sets up global state and allocates region-level data structures needed
 * for stamp building. Call this once before processing regions.
 *
 * Args:
 *   template: template image (float*, ny×nx)
 *   science: science image (float*, ny×nx)
 *   ny, nx: image dimensions
 *   n_regions_x, n_regions_y: number of regions per axis
 *   stamps_per_region_x, stamps_per_region_y: stamps per region per axis
 *
 * Returns:
 *   0 on success, -1 on allocation failure
 */
int initBuildStampsContext(float* template, float* science,
                           int ny, int nx,
                           int n_regions_x, int n_regions_y,
                           int stamps_per_region_x, int stamps_per_region_y);

/*
 * Build stamps for a single region with proper global state management.
 *
 * Handles region data extraction, mask setup, and global state configuration
 * for a single region. Returns pointers to extracted region buffers that can
 * be passed directly to the C buildStamps function.
 *
 * Args:
 *   region_x, region_y: region coordinates (0-indexed)
 *   out_template_region: output pointer to extracted template region buffer
 *   out_science_region: output pointer to extracted science region buffer
 *
 * Returns:
 *   0 on success, -1 on error
 *
 * Note: After calling this function, rPixX and rPixY globals are set to the
 * region dimensions, and mRData points to the mask array. Use these values
 * when calling buildStamps with the returned region buffers.
 */
int buildStampsRegion(int region_x, int region_y,
                      float** out_template_region, float** out_science_region);

/*
 * Clean up buildStamps context and free allocated memory.
 *
 * Call this after processing all regions to free temporary buffers
 * and restore global state.
 */
void cleanupBuildStampsContext(void);

/*
 * Set up globals for a full-image spatial_convolve() call.
 *
 * Allocates a zero-filled pixel mask of size nx*ny, sets the TLS globals
 * rPixX, rPixY, and mRData so that spatial_convolve() works correctly on
 * the full image.
 *
 * Returns:
 *   Pointer to the mask array (pass as cMask to spatial_convolve), or NULL on error.
 */
int* setupSpatialConvolve(int nx, int ny);

/*
 * Clean up after a spatial_convolve() call.
 * Resets mRData to NULL.
 */
void cleanupSpatialConvolve(void);

/*
 * Compute kernel norm (integral) from fitted coefficients at the image center.
 *
 * Sets rPixX/rPixY to nx/ny so the spatial polynomial normalises correctly,
 * then evaluates make_kernel() at the image centre (nx/2, ny/2).
 * Call this after fitKernel() has populated kernel_coeffs.
 *
 * Args:
 *   kernel_coeffs: fitted kernel coefficients from fitKernel()
 *   nx: full image width (pixels)
 *   ny: full image height (pixels)
 *
 * Returns:
 *   kernel integral (sum of all kernel pixel values at image centre)
 */
double computeKernelNorm(double* kernel_coeffs, int nx, int ny);

/* =====================================================================
 * Threading Control Functions
 * =====================================================================
 * These functions enable and query multi-threaded execution in FFTW3 and BLAS.
 * Must be called once before any kernel fitting or convolution operations.
 */

/*
 * Initialize FFTW3 and BLAS/LAPACK threading for multi-threaded execution.
 *
 * This function must be called once before any HOTPANTS kernels are executed.
 * It enables:
 *   - FFTW3 multi-threaded FFT execution (if libfftw3_omp is available)
 *   - BLAS/LAPACK multi-threaded linear algebra (OpenBLAS, MKL, etc.)
 *
 * Sets threading to the OpenMP max threads if available, otherwise 1.
 *
 * Safe to call multiple times (idempotent).
 *
 * Returns:
 *   0 on success, -1 if FFTW3 threading initialization failed
 */
int init_threading(void);

/*
 * Query the number of threads being used by BLAS/LAPACK.
 *
 * Returns:
 *   Number of threads (1 if single-threaded BLAS or unknown implementation)
 */
int get_blas_threads(void);

/*
 * Query FFTW3 threading availability.
 *
 * Returns:
 *   1 if FFTW3 threading is initialized and available, 0 otherwise
 */
int get_fftw3_threading_available(void);

/*
 * Adjust BLAS thread count based on region layout to avoid oversubscription.
 *
 * In multi-region mode, sets BLAS threads to 1 (region-level parallelism is sufficient).
 * In single-region mode, sets BLAS threads to omp_get_max_threads() (maximize per-region parallelism).
 *
 * Should be called once at the start of processing, before fitKernel().
 *
 * Args:
 *   nRegX, nRegY: Number of regions per axis
 *
 * Returns:
 *   0 on success
 */
int adjust_blas_threads_for_region_layout(int nRegX, int nRegY);

/* =====================================================================
 * Kernel Basis Type Configuration
 * =====================================================================
 * These functions manage the kernel basis type selection (Gaussian vs. Delta).
 * Call before fitKernel() or spatial_convolve() to switch basis implementations.
 */

/*
 * Set the kernel basis type.
 *
 * Args:
 *   basis_type: BASIS_TYPE_GAUSSIAN (0) or BASIS_TYPE_DELTA (1)
 *
 * Returns:
 *   0 on success, -1 if basis_type is invalid
 */
int set_basis_type(int basis_type);

/*
 * Get the current kernel basis type.
 *
 * Returns:
 *   Current basis type (BASIS_TYPE_GAUSSIAN or BASIS_TYPE_DELTA)
 */
int get_basis_type(void);

/*
 * Set delta basis grid spacing (pixels).
 *
 * Only used when basis_type == BASIS_TYPE_DELTA.
 *
 * Args:
 *   grid_size: Spacing between delta basis function centers (> 0)
 *
 * Returns:
 *   0 on success, -1 if grid_size <= 0
 */
int set_delta_ker_grid_size(double grid_size);

/*
 * Get delta basis grid spacing.
 *
 * Returns:
 *   Grid spacing in pixels
 */
double get_delta_ker_grid_size(void);

/*
 * Set delta basis Laplacian regularization weight.
 *
 * Only used when basis_type == BASIS_TYPE_DELTA.
 * Higher values enforce smoother kernels; 0 = no regularization.
 *
 * Args:
 *   regularization_weight: Regularization strength (>= 0)
 *
 * Returns:
 *   0 on success, -1 if regularization_weight < 0
 */
int set_delta_regularization(double regularization_weight);

/*
 * Get delta basis Laplacian regularization weight.
 *
 * Returns:
 *   Current regularization weight
 */
double get_delta_regularization(void);

/* =====================================================================
 * Thread-Safe _ctx Wrapper Functions
 * =====================================================================
 * Each function acquires an internal mutex, applies cfg to globals,
 * calls the core algorithm, reads back computed outputs, releases mutex.
 *
 * Preferred call sequence for fit_kernel:
 *   initKernelGlobals_ctx(cfg, nx, ny, ...)
 *   → allocateStamps()
 *   → initBuildStampsContext_ctx(cfg, ...)
 *   → for each region: buildStampsRegion_ctx + fillStampsForRegion_ctx
 *   → fitKernel_ctx(cfg, ...)
 *   → computeKernelNorm_ctx(cfg, ...)
 *   → cleanupBuildStampsContext_ctx()
 *
 * Preferred call sequence for spatial_convolve:
 *   initKernelGlobals_ctx(cfg, nx, ny, ...)
 *   → setupSpatialConvolve_ctx(cfg, nx, ny)
 *   → spatial_convolve_ctx(cfg, ...)
 *   → cleanupSpatialConvolve_ctx()
 */

int    initKernelGlobals_ctx(hotpants_config_t* cfg,
                              int image_nx, int image_ny,
                              int n_reg_x, int n_reg_y,
                              int n_stamp_x, int n_stamp_y);

int    initBuildStampsContext_ctx(hotpants_config_t* cfg,
                                   float* tmpl, float* science,
                                   int ny, int nx,
                                   int n_regions_x, int n_regions_y,
                                   int stamps_per_x, int stamps_per_y);

int    buildStampsRegion_ctx(hotpants_config_t* cfg,
                              int region_x, int region_y,
                              float** out_template_region,
                              float** out_science_region);

int    fillStampsForRegion_ctx(hotpants_config_t* cfg,
                                stamp_struct* stamps, int n_stamps);

void   fitKernel_ctx(hotpants_config_t* cfg,
                     stamp_struct* stamps,
                     float* imRef, float* imConv, float* imNoise,
                     double* kernel_coeffs,
                     double* meansig, double* scatter, int* n_skipped);

int*   setupSpatialConvolve_ctx(hotpants_config_t* cfg, int nx, int ny);

double computeKernelNorm_ctx(hotpants_config_t* cfg,
                              double* kernel_coeffs, int nx, int ny);

void   spatial_convolve_ctx(hotpants_config_t* cfg,
                             float* image, float** var_image,
                             int ny, int nx,
                             double* kernel_coeffs, float* output,
                             int* conv_method);

void   cleanupBuildStampsContext_ctx(void);
void   cleanupSpatialConvolve_ctx(void);

#endif /* HOTPANTS_API_H */
