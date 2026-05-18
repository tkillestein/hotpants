#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "defaults.h"
#include "globals.h"
#include "functions.h"

/* BLAS/LAPACK threading control macros */
#ifdef __OPENBLAS__
extern int openblas_get_num_threads(void);
extern void openblas_set_num_threads(int num_threads);
#endif

#ifdef __MKL__
#include <mkl.h>
#endif

/* Track whether FFTW3 threading has been initialized (thread-safe: called once) */
static int fftw3_threading_initialized = 0;

/**
 * @brief Initialize FFTW3 and BLAS threading for multi-threaded execution.
 *
 * This function must be called once before any HOTPANTS kernels are executed.
 * It enables:
 *   - FFTW3 multi-threaded FFT execution (if fftw3_omp is available)
 *   - BLAS/LAPACK multi-threaded linear algebra (OpenBLAS or MKL)
 *
 * Sets threading to the OpenMP max threads if available, otherwise 1.
 *
 * Safe to call multiple times (idempotent).
 *
 * @return 0 on success, -1 if unable to initialize FFTW3 threads
 */
int init_threading(void) {
  if (fftw3_threading_initialized)
    return 0;  /* Already initialized */

  int num_threads = 1;
#ifdef _OPENMP
  num_threads = omp_get_max_threads();
#endif

  /* Initialize FFTW3 multi-threaded execution.
   * fftw_init_threads() prepares the library for threaded operation.
   * fftw_plan_with_nthreads(n) tells FFTW3 to use n threads for plan creation and execution.
   * Must be called BEFORE any fftw_plan_*() calls.
   */
  if (fftw_init_threads() == 0) {
    fprintf(stderr, "WARNING: fftw_init_threads() failed; FFT will be single-threaded\n");
    return -1;
  }
  fftw_plan_with_nthreads(num_threads);

  /* Initialize BLAS/LAPACK threading.
   * Different BLAS implementations have different APIs; we support OpenBLAS and MKL.
   */
#ifdef __OPENBLAS__
  openblas_set_num_threads(num_threads);
#endif

#ifdef __MKL__
  mkl_set_num_threads(num_threads);
#endif

  fftw3_threading_initialized = 1;
  return 0;
}

/**
 * @brief Query the number of threads being used by BLAS/LAPACK.
 *
 * @return Number of threads (1 if single-threaded BLAS or unknown implementation)
 */
int get_blas_threads(void) {
#ifdef __OPENBLAS__
  return openblas_get_num_threads();
#elif defined(__MKL__)
  return mkl_get_max_threads();
#else
  return 1;  /* Single-threaded or unknown */
#endif
}

/**
 * @brief Get FFTW3 threading status.
 *
 * @return 1 if FFTW3 threading is initialized and available, 0 otherwise
 */
int get_fftw3_threading_available(void) {
  return fftw3_threading_initialized;
}

/**
 * @brief Adjust BLAS thread count based on region layout to avoid oversubscription.
 *
 * In multi-region mode (nRegX > 1 or nRegY > 1), region-level OpenMP parallelism
 * provides sufficient parallelism, so BLAS is set to single-threaded to avoid spawning
 * too many threads (e.g., 4 regions × 4 BLAS threads = 16 threads on 4-core system).
 *
 * In single-region mode, BLAS is configured to use all available threads for maximum
 * parallelism within the Cholesky solve and other linear algebra operations.
 *
 * @param nRegX, nRegY: Number of regions per axis
 * @return 0 on success
 */
int adjust_blas_threads_for_region_layout(int nRegX, int nRegY) {
  int blas_threads = 1;

  /* Single-region mode: use all threads for BLAS */
  if (nRegX == 1 && nRegY == 1) {
#ifdef _OPENMP
    blas_threads = omp_get_max_threads();
#else
    blas_threads = 1;
#endif
  }
  /* Multi-region mode: use single thread for BLAS (region loop provides parallelism) */

#ifdef __OPENBLAS__
  openblas_set_num_threads(blas_threads);
#endif

#ifdef __MKL__
  mkl_set_num_threads(blas_threads);
#endif

  return 0;
}

/* Declare globals accessed by wrapper */
extern __thread int rPixX, rPixY;
extern __thread int* mRData;
extern char* forceConvolve, *photNormalize, *figMerit;
extern int fwKernel, fwStamp;
extern int fwKSStamp;
extern int kcStep;

/* Basis globals (set in vargs.c but needed before allocateStamps) */
extern int ngauss, *deg_fixe;
extern float *sigma_gauss;
extern int nCompKer, nComp, nC, nCompBG, nBGVectors, nCompTotal;
extern int hwKernel, kerOrder, bgOrder;
extern int nKSStamps, hwKSStamp;
extern int nStampX, nStampY;
extern int nRegX, nRegY;

/* Algorithm tuning globals normally set by vargs() */
extern float statSig, kerSigReject, kerFracMask;
extern float kfSpreadMask1, kfSpreadMask2;
extern float fillVal, fillValNoise;

/* Kernel basis arrays — allocated once per parameter set, used by fillStamp
 * and spatial_convolve (mirrors main.c lines 450-461). */
extern double **kernel_vec;
extern double *filter_x, *filter_y;
/* TLS kernel scratch buffers for convolution (main thread copies) */
extern __thread double *kernel, *kernel_coeffs;
/* TLS scratch buffer used by make_model / getStampSig inside fitKernel */
extern __thread float *temp;

/* Forward declarations - defined in alard.c */
extern int getKernelVec(void);
extern double make_kernel(int xi, int yi, double* kernelSol);
extern double make_kernel_dispatch(int xi, int yi, double* kernelSol);

/* Forward declarations - defined in functions.c */
extern void makeInputMask(float* iRData, float* tRData, int* mRData);

/* Forward declaration - defined in alard.c */
extern int fillStamp(stamp_struct* stamp, float* imConv, float* imRef);

/* Scratch mask buffer for spatial_convolve when no region context is active */
static struct {
  int* mask_buffer;
  int mask_size;
} convolve_context = { NULL, 0 };

/* Wrapper context: global state for buildStamps operations */
static struct {
  float* template_image;      /* Full template image (read-only) */
  float* science_image;       /* Full science image (read-only) */
  int full_ny, full_nx;       /* Full image dimensions */
  int n_regions_x, n_regions_y;
  int stamps_per_region_x, stamps_per_region_y;

  float* region_t_buffer;     /* Per-region template buffer */
  float* region_i_buffer;     /* Per-region science buffer */
  int* region_mask_buffer;    /* Per-region mask array */

  int initialized;            /* Flag to track initialization */
} wrapper_context = {
  .initialized = 0
};

/*
 * Initialize kernel basis globals from hwKernel/kerOrder/bgOrder.
 *
 * This replicates the setup done in main.c lines 400-441 so that allocateStamps
 * and the C fitting code have correctly-sized arrays. Must be called before
 * allocateStamps() and buildStamps().
 *
 * Uses default 3-Gaussian basis (matching D_NGAUSS, D_DEG_GAUSSn, D_SIG_GAUSSn
 * from defaults.h):
 *   Gaussian 0: deg=6, sigma=0.7 px
 *   Gaussian 1: deg=4, sigma=1.5 px
 *   Gaussian 2: deg=2, sigma=3.0 px
 *
 * Args:
 *   image_nx, image_ny: full image dimensions (used to compute fwStamp)
 *   n_reg_x, n_reg_y: number of regions (used to compute fwStamp)
 *   n_stamp_x, n_stamp_y: stamps per region (used to compute fwStamp)
 *
 * Returns:
 *   0 on success, -1 on allocation failure
 */
int initKernelGlobals(int image_nx, int image_ny,
                      int n_reg_x, int n_reg_y,
                      int n_stamp_x, int n_stamp_y) {
  int i;
  /* Tracks how many entries were allocated in kernel_vec so cleanup is safe
   * across calls where iBasisType may differ (delta overrides nCompKer). */
  static int kernel_vec_alloc_size = 0;

  /* --- Gaussian basis setup (mirrors vargs.c lines 36-45) --- */

  /* Free any existing allocations */
  if (deg_fixe)   { free(deg_fixe);   deg_fixe = NULL; }
  if (sigma_gauss){ free(sigma_gauss); sigma_gauss = NULL; }

  ngauss = 3; /* D_NGAUSS */
  deg_fixe = (int*)calloc(ngauss, sizeof(int));
  sigma_gauss = (float*)calloc(ngauss, sizeof(float));
  if (!deg_fixe || !sigma_gauss) return -1;

  /* Default Gaussian degrees (D_DEG_GAUSS1/2/3) */
  deg_fixe[0] = 6;
  deg_fixe[1] = 4;
  deg_fixe[2] = 2;

  /* sigma_gauss stored as 1/(2*sigma^2) — same convention as vargs.c */
  sigma_gauss[0] = 1.0f / (2.0f * 0.7f * 0.7f);   /* D_SIG_GAUSS1 = 0.7 */
  sigma_gauss[1] = 1.0f / (2.0f * 1.5f * 1.5f);   /* D_SIG_GAUSS2 = 1.5 */
  sigma_gauss[2] = 1.0f / (2.0f * 3.0f * 3.0f);   /* D_SIG_GAUSS3 = 3.0 */

  /* --- Derived counts for Gaussian basis (mirrors main.c lines 401-409) --- */

  nCompKer = 0;
  for (i = 0; i < ngauss; i++)
    nCompKer += ((deg_fixe[i] + 1) * (deg_fixe[i] + 2)) / 2;

  nComp   = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  nBGVectors = ((bgOrder + 1) * (bgOrder + 2)) / 2;

  /* --- Frame widths (mirrors main.c lines 411-440) --- */

  fwKernel  = hwKernel * 2 + 1;
  fwKSStamp = hwKSStamp * 2 + 1;

  /* fwStamp: stamp frame width derived from image + region + stamp grid */
  int min_dim_per_region_x = image_nx / n_reg_x / n_stamp_x;
  int min_dim_per_region_y = image_ny / n_reg_y / n_stamp_y;
  int min_dim = (min_dim_per_region_x < min_dim_per_region_y)
                ? min_dim_per_region_x : min_dim_per_region_y;
  fwStamp = min_dim - fwKernel;
  if (fwStamp % 2 == 0) fwStamp -= 1;  /* keep odd */

  /* Sanity-check: stamp must be at least as large as substamp */
  if (fwStamp < fwKSStamp) {
    fwStamp = fwKSStamp + fwKernel;
    if (fwStamp % 2 == 0) fwStamp -= 1;
    fprintf(stderr, "WARNING: fwStamp < fwKSStamp, adjusted to %d\n", fwStamp);
  }

  /* kcStep: convolution step size (mirrors main.c line 442) */
  kcStep = fwKernel;

  /* --- Kernel basis array allocation (mirrors main.c lines 450-461) ---
   * kernel_vec, filter_x, filter_y hold the pre-computed Gaussian-polynomial
   * basis images and separable filters used by fillStamp/xy_conv_stamp and
   * spatial_convolve/make_kernel_local.
   * Use kernel_vec_alloc_size (not nCompKer) for cleanup to handle re-entry
   * after a delta basis run that overrides nCompKer below. */

  if (kernel_vec) {
    for (i = 0; i < kernel_vec_alloc_size; i++) {
      if (kernel_vec[i]) { free(kernel_vec[i]); kernel_vec[i] = NULL; }
    }
    free(kernel_vec); kernel_vec = NULL;
  }
  if (filter_x) { free(filter_x); filter_x = NULL; }
  if (filter_y) { free(filter_y); filter_y = NULL; }
  if (kernel)        { free(kernel);        kernel = NULL; }
  if (kernel_coeffs) { free(kernel_coeffs); kernel_coeffs = NULL; }
  if (temp)          { free(temp);          temp = NULL; }

  /* Allocate using Gaussian nCompKer; delta basis does not use these arrays. */
  kernel_vec_alloc_size = nCompKer;
  kernel_vec = (double**)calloc(nCompKer, sizeof(double*));
  filter_x   = (double*)calloc(nCompKer * fwKernel, sizeof(double));
  filter_y   = (double*)calloc(nCompKer * fwKernel, sizeof(double));
  kernel     = (double*)calloc(fwKernel * fwKernel, sizeof(double));
  kernel_coeffs = (double*)calloc(nCompKer, sizeof(double));
  /* temp: scratch buffer for make_model/getStampSig (mirrors main.c line 450) */
  temp = (float*)calloc((fwKSStamp + fwKernel) * fwKSStamp, sizeof(float));

  if (!kernel_vec || !filter_x || !filter_y || !kernel || !kernel_coeffs || !temp)
    return -1;

  /* Populate kernel_vec[i] and filter_x/filter_y via getKernelVec().
   * This replicates the per-region call at main.c lines 1093-1094. */
  int nbasis = getKernelVec();
  if (nbasis < 0) {
    fprintf(stderr, "ERROR: Failed to initialize Gaussian kernel basis\n");
    return -1;
  }

  /* --- Delta basis: override nCompKer = fwKernel² ---
   * For delta basis, each kernel pixel is a separate basis function, so
   * nCompKer must equal fwKernel² (not the Gaussian polynomial count).
   * This override must happen BEFORE allocateStamps() is called so that
   * stamp->mat, stamp->scprod, and stamp->vectors are sized correctly.
   * kernel_vec/filter_x/filter_y are Gaussian-only; delta fillStamp uses
   * xy_conv_stamp_delta() which does not touch these arrays. */
  if (iBasisType == BASIS_TYPE_DELTA) {
    nCompKer = fwKernel * fwKernel;
    fprintf(stderr, "INFO: Delta basis active: nCompKer = %d (%d×%d kernel pixels)\n",
            nCompKer, fwKernel, fwKernel);
  }

  /* Recompute derived counts with the (possibly overridden) nCompKer */
  nC      = nCompKer + 2;  /* used by allocateStamps for mat/scprod size */
  nCompBG = (nCompKer - 1) * nComp + 1;
  nCompTotal = nCompKer * nComp + nBGVectors;

  /* Initialize string globals normally set by vargs() — these must not be NULL
   * because buildStamps, fitKernel, and spatial_convolve call strncmp() on them.
   * forceConvolve = "t": Python API always fits template ⊗ K ≈ science,
   * so buildStamps must only process template stamps (not both directions). */
  figMerit      = "v";   /* D_FIGMERIT: variance-based figure of merit */
  photNormalize = "t";   /* D_NORMALIZE: output on template photometric system */
  forceConvolve = "t";   /* Python API always convolves template; mirrors CLI -c t */

  /* Initialize algorithm-tuning globals normally set by vargs(); without this
   * statSig=0 and kerSigReject=0, which causes sigma_clip to clip everything
   * and check_again to reject all stamps on every iteration. */
  if (statSig == 0.0f)        statSig        = D_STATSIG;
  if (kerSigReject == 0.0f)   kerSigReject   = D_KSIGREJECT;
  if (fillVal == 0.0f)        fillVal        = D_FILL;
  if (kfSpreadMask1 == 0.0f)  kfSpreadMask1  = D_INMASKFSPREAD;
  if (kfSpreadMask2 == 0.0f)  kfSpreadMask2  = D_OUMASKFSPREAD;

  /* Set the active basis so make_kernel_dispatch() routes correctly.
   * For Gaussian this re-runs gaussian_init() (harmless since kernel_vec is
   * already populated above); for delta it runs delta_init() which sets
   * nCompKer = fwKernel² (consistent with the override above). */
  if (initKernelBasis() < 0) {
    fprintf(stderr, "ERROR: Failed to initialize kernel basis (iBasisType=%d)\n",
            iBasisType);
    return -1;
  }

  return 0;
}

/* Helper: compute region boundaries */
static void get_region_bounds(int region_x, int region_y, int n_regions_x,
                              int n_regions_y, int full_nx, int full_ny,
                              int* rx_min, int* rx_max, int* ry_min, int* ry_max,
                              int* r_pix_x, int* r_pix_y) {
  /* Evenly divide image into regions */
  *rx_min = (region_x * full_nx) / n_regions_x;
  *rx_max = ((region_x + 1) * full_nx) / n_regions_x - 1;
  *ry_min = (region_y * full_ny) / n_regions_y;
  *ry_max = ((region_y + 1) * full_ny) / n_regions_y - 1;

  *r_pix_x = *rx_max - *rx_min + 1;
  *r_pix_y = *ry_max - *ry_min + 1;
}

/* Initialize buildStamps context */
int initBuildStampsContext(float* template, float* science,
                           int ny, int nx,
                           int n_regions_x, int n_regions_y,
                           int stamps_per_region_x, int stamps_per_region_y) {
  if (wrapper_context.initialized) {
    fprintf(stderr, "WARNING: Context already initialized\n");
    return -1;
  }

  /* Validate inputs */
  if (!template || !science || ny <= 0 || nx <= 0) {
    fprintf(stderr, "ERROR: Invalid inputs to initBuildStampsContext\n");
    return -1;
  }

  /* Save image pointers and dimensions */
  wrapper_context.template_image = template;
  wrapper_context.science_image = science;
  wrapper_context.full_ny = ny;
  wrapper_context.full_nx = nx;
  wrapper_context.n_regions_x = n_regions_x;
  wrapper_context.n_regions_y = n_regions_y;
  wrapper_context.stamps_per_region_x = stamps_per_region_x;
  wrapper_context.stamps_per_region_y = stamps_per_region_y;

  /* Estimate maximum region size for buffer allocation.
   * Add hwKernel border on each side to match main.c region extraction. */
  int max_region_x = (nx + n_regions_x - 1) / n_regions_x + 1 + 2 * hwKernel;
  int max_region_y = (ny + n_regions_y - 1) / n_regions_y + 1 + 2 * hwKernel;
  int max_region_size = max_region_x * max_region_y;

  /* Allocate region buffers */
  wrapper_context.region_t_buffer = (float*)calloc(max_region_size, sizeof(float));
  wrapper_context.region_i_buffer = (float*)calloc(max_region_size, sizeof(float));
  wrapper_context.region_mask_buffer = (int*)calloc(max_region_size, sizeof(int));

  if (!wrapper_context.region_t_buffer || !wrapper_context.region_i_buffer ||
      !wrapper_context.region_mask_buffer) {
    fprintf(stderr, "ERROR: Failed to allocate region buffers\n");
    return -1;
  }

  /* Set global forceConvolve to "t": fit template ⊗ K ≈ science.
   * This is the natural direction when science has a broader PSF than template,
   * and produces a difference image D = science - template ⊗ K. */
  forceConvolve = "t";

  wrapper_context.initialized = 1;
  return 0;
}

/* Extract and setup region data for buildStamps
 *
 * This function extracts a region from the full images, sets up global state,
 * and prepares data for the C buildStamps function. Python code then calls
 * buildStamps directly with the returned pointers.
 */
int buildStampsRegion(int region_x, int region_y,
                      float** out_template_region, float** out_science_region) {
  if (!wrapper_context.initialized) {
    fprintf(stderr, "ERROR: Context not initialized\n");
    return -1;
  }

  /* Compute region boundaries (without border) */
  int rx_min, rx_max, ry_min, ry_max;
  int r_pix_x, r_pix_y;
  get_region_bounds(region_x, region_y,
                    wrapper_context.n_regions_x, wrapper_context.n_regions_y,
                    wrapper_context.full_nx, wrapper_context.full_ny,
                    &rx_min, &rx_max, &ry_min, &ry_max,
                    &r_pix_x, &r_pix_y);

  /* Add hwKernel border, clamped to image bounds — mirrors main.c lines 709-723.
   * buildStamps/fillStamp require rPixX to be the bordered width and iRData to
   * start at rXBMin = rx_min - hwKernel; without this the buffer is indexed
   * out-of-bounds and the process segfaults. */
  int hw = hwKernel;
  int rx_b_min = (rx_min >= hw) ? rx_min - hw : 0;
  int rx_b_max = (rx_max + hw < wrapper_context.full_nx)
                   ? rx_max + hw : wrapper_context.full_nx - 1;
  int ry_b_min = (ry_min >= hw) ? ry_min - hw : 0;
  int ry_b_max = (ry_max + hw < wrapper_context.full_ny)
                   ? ry_max + hw : wrapper_context.full_ny - 1;
  int r_b_pix_x = rx_b_max - rx_b_min + 1;
  int r_b_pix_y = ry_b_max - ry_b_min + 1;

  /* Set global region dimensions to the bordered size (required by buildStamps/fillStamp) */
  rPixX = r_b_pix_x;
  rPixY = r_b_pix_y;

  /* Extract region data WITH border from full images */
  for (int y = 0; y < r_b_pix_y; y++) {
    for (int x = 0; x < r_b_pix_x; x++) {
      int region_idx = y * r_b_pix_x + x;
      int full_idx = (ry_b_min + y) * wrapper_context.full_nx + (rx_b_min + x);
      wrapper_context.region_t_buffer[region_idx] = wrapper_context.template_image[full_idx];
      wrapper_context.region_i_buffer[region_idx] = wrapper_context.science_image[full_idx];
    }
  }

  /* Initialize mask array for bordered region */
  memset(wrapper_context.region_mask_buffer, 0, r_b_pix_x * r_b_pix_y * sizeof(int));

  /* Set the global mask array pointer */
  mRData = wrapper_context.region_mask_buffer;

  /* Set up mask from image data (template first, then science — matches main.c) */
  makeInputMask(wrapper_context.region_t_buffer, wrapper_context.region_i_buffer,
                wrapper_context.region_mask_buffer);

  /* Apply border masking — mirrors main.c lines 893-907.
   * sBorder = hwKernel + kcStep pixels at each edge are flagged FLAG_T_BAD|FLAG_I_BAD.
   * This ensures getPsfCenters never places a PSF center within sBorder of the
   * buffer edge, guaranteeing xy_conv_stamp won't access out-of-bounds memory
   * (which requires hwKSStamp+hwKernel pixels of margin from any PSF center). */
  int sBorder = hwKernel + kcStep;
  for (int by = 0; by < r_b_pix_y; by++) {
    for (int bx = 0; bx < sBorder; bx++)
      wrapper_context.region_mask_buffer[bx + r_b_pix_x * by] |= (FLAG_T_BAD | FLAG_I_BAD);
    for (int bx = r_b_pix_x - sBorder; bx < r_b_pix_x; bx++)
      wrapper_context.region_mask_buffer[bx + r_b_pix_x * by] |= (FLAG_T_BAD | FLAG_I_BAD);
  }
  for (int by = 0; by < sBorder; by++)
    for (int bx = sBorder; bx < r_b_pix_x - sBorder; bx++)
      wrapper_context.region_mask_buffer[bx + r_b_pix_x * by] |= (FLAG_T_BAD | FLAG_I_BAD);
  for (int by = r_b_pix_y - sBorder; by < r_b_pix_y; by++)
    for (int bx = sBorder; bx < r_b_pix_x - sBorder; bx++)
      wrapper_context.region_mask_buffer[bx + r_b_pix_x * by] |= (FLAG_T_BAD | FLAG_I_BAD);

  /* Return pointers to extracted region data */
  *out_template_region = wrapper_context.region_t_buffer;
  *out_science_region = wrapper_context.region_i_buffer;

  return 0; /* Success */
}

/*
 * Set up globals for spatial_convolve on a full image.
 *
 * spatial_convolve() uses:
 *   rPixX, rPixY  — image dimensions for kernel polynomial normalisation
 *   mRData        — pixel mask array (updated with convolution flags)
 *
 * This function allocates a zero-filled mask of size nx*ny, sets mRData and
 * rPixX/rPixY, and returns a pointer to the mask (passed as cMask to
 * spatial_convolve).
 *
 * Args:
 *   nx, ny: full image dimensions (width, height)
 *
 * Returns:
 *   Pointer to allocated int mask array on success, NULL on failure.
 *   The pointer is valid until cleanupSpatialConvolve() is called.
 */
int* setupSpatialConvolve(int nx, int ny) {
  int new_size = nx * ny;
  if (new_size <= 0) return NULL;

  if (convolve_context.mask_buffer && convolve_context.mask_size != new_size) {
    free(convolve_context.mask_buffer);
    convolve_context.mask_buffer = NULL;
  }

  if (!convolve_context.mask_buffer) {
    convolve_context.mask_buffer = (int*)calloc(new_size, sizeof(int));
    if (!convolve_context.mask_buffer) return NULL;
    convolve_context.mask_size = new_size;
  } else {
    memset(convolve_context.mask_buffer, 0, new_size * sizeof(int));
  }

  /* Broadcast TLS variables to all OMP threads.  spatial_convolve runs an
   * #pragma omp parallel block, so each worker thread needs its own copy of
   * mRData, rPixX, rPixY set before entering that region.  Setting them only
   * in the main thread would leave worker-thread copies NULL/0. */
  int* buf = convolve_context.mask_buffer;
#ifdef _OPENMP
#pragma omp parallel
  {
    mRData = buf;
    rPixX  = nx;
    rPixY  = ny;
  }
#else
  mRData = buf;
  rPixX  = nx;
  rPixY  = ny;
#endif

  return convolve_context.mask_buffer;
}

/*
 * Clean up after spatial_convolve call.
 * Resets mRData to NULL in all OMP threads.
 */
void cleanupSpatialConvolve(void) {
#ifdef _OPENMP
#pragma omp parallel
  {
    mRData = NULL;
  }
#else
  mRData = NULL;
#endif
}

/*
 * Fill stamp vectors for kernel fitting.
 *
 * Must be called after buildStamps() for a region, and before fitKernel(),
 * while the region context (mRData, region buffers, rPixX/rPixY) is still
 * active.  Mirrors main.c lines 1099-1119.
 *
 * Args:
 *   stamps: stamp array pointer (element 0 is the first stamp for this region)
 *   n_stamps: number of stamps to fill
 *
 * Returns 0 on success, -1 if context is not initialized.
 */
int fillStampsForRegion(stamp_struct* stamps, int n_stamps) {
  if (!wrapper_context.initialized) {
    LOG_ERROR("fillStampsForRegion: context not initialized");
    return -1;
  }
  int k;
  for (k = 0; k < n_stamps; k++) {
    /* forceConvolve = "t": imConv = template, imRef = science.
     * Fits template ⊗ K ≈ science; difference = science - template ⊗ K. */
    fillStamp(&stamps[k],
              wrapper_context.region_t_buffer,
              wrapper_context.region_i_buffer);
  }
  return 0;
}

/* Clean up buildStamps context */
void cleanupBuildStampsContext(void) {
  if (!wrapper_context.initialized) {
    return;
  }

  /* Free allocated buffers */
  if (wrapper_context.region_t_buffer) {
    free(wrapper_context.region_t_buffer);
    wrapper_context.region_t_buffer = NULL;
  }
  if (wrapper_context.region_i_buffer) {
    free(wrapper_context.region_i_buffer);
    wrapper_context.region_i_buffer = NULL;
  }
  if (wrapper_context.region_mask_buffer) {
    free(wrapper_context.region_mask_buffer);
    wrapper_context.region_mask_buffer = NULL;
  }

  /* Reset global pointers */
  mRData = NULL;

  /* Mark context as uninitialized */
  wrapper_context.initialized = 0;
}

/*
 * Evaluate kernel norm (integral) from fitted coefficients at image center.
 *
 * Sets rPixX/rPixY to full image dimensions so the spatial polynomial
 * evaluates correctly, then calls make_kernel at the image center.
 * Call this after fitKernel() has populated kernel_coeffs.
 *
 * Args:
 *   kernel_coeffs: fitted coefficients from fitKernel()
 *   nx: full image width (pixels)
 *   ny: full image height (pixels)
 *
 * Returns:
 *   kernel integral (sum of all kernel pixel values at image center)
 */
double computeKernelNorm(double* kernel_coeffs, int nx, int ny) {
  rPixX = nx;
  rPixY = ny;
  /* Use the dispatch function so delta basis routes to delta_eval_kernel_dispatch()
   * rather than the Gaussian-specific make_kernel() which accesses kernel_vec. */
  return make_kernel_dispatch(nx / 2, ny / 2, kernel_coeffs);
}

/* =====================================================================
 * Kernel Basis Type Configuration
 * =====================================================================
 * Functions for managing kernel basis type selection.
 */

/**
 * @brief Set the kernel basis type.
 *
 * @param[in] basis_type BASIS_TYPE_GAUSSIAN (0) or BASIS_TYPE_DELTA (1)
 * @return 0 on success, -1 if basis_type is invalid
 */
int set_basis_type(int basis_type) {
  if (basis_type != BASIS_TYPE_GAUSSIAN && basis_type != BASIS_TYPE_DELTA) {
    fprintf(stderr, "ERROR: Invalid basis_type %d (must be %d or %d)\n",
            basis_type, BASIS_TYPE_GAUSSIAN, BASIS_TYPE_DELTA);
    return -1;
  }
  iBasisType = basis_type;
  return 0;
}

/**
 * @brief Get the current kernel basis type.
 *
 * @return Current basis type (BASIS_TYPE_GAUSSIAN or BASIS_TYPE_DELTA)
 */
int get_basis_type(void) {
  return iBasisType;
}

/**
/**
 * @brief Set delta basis Laplacian regularization weight.
 *
 * Only used when basis_type == BASIS_TYPE_DELTA.
 * Higher values enforce smoother kernels; 0 = no regularization.
 *
 * @param[in] regularization_weight Regularization strength (>= 0)
 * @return 0 on success, -1 if regularization_weight < 0
 */
int set_delta_regularization(double regularization_weight) {
  if (regularization_weight < 0.0) {
    fprintf(stderr,
            "ERROR: deltaRegularization must be non-negative, got %f\n",
            regularization_weight);
    return -1;
  }
  rDeltaRegularization = regularization_weight;
  return 0;
}

/**
 * @brief Get delta basis Laplacian regularization weight.
 *
 * @return Current regularization weight
 */
double get_delta_regularization(void) {
  return rDeltaRegularization;
}
