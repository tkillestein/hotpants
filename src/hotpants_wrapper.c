#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "defaults.h"
#include "hotpants_context.h"

/* stamp_struct definition (copied from globals.h) */
typedef struct {
  int x0, y0;
  int x, y;
  int nx, ny;
  int* xss;
  int* yss;
  int nss;
  int sscnt;
  double** vectors;
  double* krefArea;
  double** mat;
  double* scprod;
  double sum;
  double mean;
  double median;
  double mode;
  double sd;
  double fwhm;
  double lfwhm;
  double chi2;
  double norm;
  double diff;
} stamp_struct;

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
extern void getKernelVec(void);
extern double make_kernel(int xi, int yi, double* kernelSol);
extern void fitKernel(stamp_struct* stamps, float* imRef, float* imConv,
                      float* imNoise, double* kernel_coeffs, double* meansig,
                      double* scatter, int* n_skipped);
extern void spatial_convolve(float* image, float** var_image, int ny, int nx,
                             double* kernel_coeffs, float* output, int* conv_method);

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

  /* --- Derived counts (mirrors main.c lines 401-409) --- */

  nCompKer = 0;
  for (i = 0; i < ngauss; i++)
    nCompKer += ((deg_fixe[i] + 1) * (deg_fixe[i] + 2)) / 2;

  nComp   = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  nC      = nCompKer + 2;  /* used by allocateStamps for mat size */
  nCompBG = (nCompKer - 1) * nComp + 1;
  nBGVectors = ((bgOrder + 1) * (bgOrder + 2)) / 2;
  nCompTotal = nCompKer * nComp + nBGVectors;

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
   * spatial_convolve/make_kernel_local. */

  if (kernel_vec) {
    for (i = 0; i < nCompKer; i++) {
      if (kernel_vec[i]) { free(kernel_vec[i]); kernel_vec[i] = NULL; }
    }
    free(kernel_vec); kernel_vec = NULL;
  }
  if (filter_x) { free(filter_x); filter_x = NULL; }
  if (filter_y) { free(filter_y); filter_y = NULL; }
  if (kernel)        { free(kernel);        kernel = NULL; }
  if (kernel_coeffs) { free(kernel_coeffs); kernel_coeffs = NULL; }
  if (temp)          { free(temp);          temp = NULL; }

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
  getKernelVec();

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
  return make_kernel(nx / 2, ny / 2, kernel_coeffs);
}

/* =====================================================================
 * Context-Aware Wrapper Functions (NEW API)
 * =====================================================================
 * These functions work with the new context structures (config_struct,
 * kernel_context_struct, region_context_struct) to eliminate global state.
 */

/*
 * Set threadprivate TLS globals from kernel context.
 *
 * Call this before invoking core algorithms that read threadprivate globals.
 * Used internally by wrapper functions.
 */
static void set_kernel_context_globals(const kernel_context_struct* kctx,
                                       const config_struct* config) {
  if (!kctx || !config) return;

  /* Set component counts (read by allocateStamps, fillStamp) */
  nCompKer = kctx->basis.n_comp_ker;
  nComp = kctx->basis.n_comp;
  nC = kctx->basis.n_c;
  nCompBG = kctx->basis.n_comp_bg;
  nBGVectors = kctx->basis.n_bg_vectors;
  nCompTotal = kctx->basis.n_comp_total;

  /* Set frame widths (read by fillStamp, spatial_convolve) */
  fwKernel = kctx->fw_kernel;
  fwStamp = kctx->fw_stamp;
  fwKSStamp = kctx->fw_ks_stamp;
  kcStep = kctx->kc_step;

  /* Set Gaussian basis (read by fillStamp) */
  ngauss = kctx->basis.ngauss;
  deg_fixe = kctx->basis.deg_fixe;
  sigma_gauss = kctx->basis.sigma_gauss;

  /* Set kernel basis vectors (read by fillStamp, spatial_convolve) */
  kernel_vec = kctx->basis.kernel_vec;
  filter_x = kctx->basis.filter_x;
  filter_y = kctx->basis.filter_y;

  /* Set algorithm tuning globals (read by fitKernel) */
  statSig = config->stat_sig;
  kerSigReject = config->ker_sig_reject;
  kerFracMask = config->ker_frac_mask;
  fillVal = config->fill_val;
  fillValNoise = config->fill_val_noise;
  kfSpreadMask1 = config->kf_spread_mask1;
  kfSpreadMask2 = config->kf_spread_mask2;

  /* Set configuration kernel parameters (read by buildStamps) */
  hwKernel = config->kernel_half_width;
  kerOrder = config->kernel_order;
  bgOrder = config->bg_order;
  nKSStamps = config->n_ks_stamps;
  hwKSStamp = config->hw_ks_stamp;

  /* Set string options (read by buildStamps, fitKernel) */
  static char force_buf[2] = {0, 0};
  force_buf[0] = config->force_convolve;
  forceConvolve = force_buf;
  figMerit = "v";
  photNormalize = "t";
}

/*
 * Set threadprivate region context globals.
 *
 * Call this before invoking core algorithms on a region.
 * Used internally by wrapper functions.
 */
static void set_region_context_globals(const region_context_struct* rctx) {
  if (!rctx) return;

  rPixX = rctx->pix_x;
  rPixY = rctx->pix_y;
  mRData = rctx->mask;
  kernel = rctx->kernel;
  kernel_coeffs = rctx->kernel_coeffs;
  temp = rctx->temp;
}

/*
 * Broadcast threadprivate context globals to all OpenMP threads.
 *
 * Use inside OpenMP parallel regions to ensure all threads have
 * correct copies of TLS variables.
 */
static void broadcast_region_context_globals(const region_context_struct* rctx) {
  if (!rctx) return;

#ifdef _OPENMP
#pragma omp parallel
  {
    set_region_context_globals(rctx);
  }
#else
  set_region_context_globals(rctx);
#endif
}

/*
 * Initialize context-aware buildStamps operation.
 *
 * Sets up kernel and region globals from contexts, preparing for
 * a buildStamps call with the given configuration.
 *
 * Args:
 *   config: configuration struct
 *   kctx: kernel context (must be initialized)
 *   rctx: region context (must be initialized)
 *
 * Returns: 0 on success, -1 on error
 */
int hotpants_setup_buildstamps(const config_struct* config,
                               const kernel_context_struct* kctx,
                               region_context_struct* rctx) {
  if (!config || !kctx || !rctx) return -1;

  set_kernel_context_globals(kctx, config);
  set_region_context_globals(rctx);

  return 0;
}

/*
 * Initialize context-aware fitKernel operation.
 *
 * Sets up kernel and region globals from contexts, preparing for
 * a fitKernel call.
 *
 * Args:
 *   config: configuration struct
 *   kctx: kernel context (must be initialized)
 *   rctx: region context (must be initialized)
 *
 * Returns: 0 on success, -1 on error
 */
int hotpants_setup_fitkernel(const config_struct* config,
                             const kernel_context_struct* kctx,
                             region_context_struct* rctx) {
  if (!config || !kctx || !rctx) return -1;

  set_kernel_context_globals(kctx, config);
  broadcast_region_context_globals(rctx);

  return 0;
}

/*
 * Initialize context-aware spatial_convolve operation.
 *
 * Sets up kernel and region globals from contexts, preparing for
 * a spatial_convolve call.
 *
 * Args:
 *   config: configuration struct
 *   kctx: kernel context (must be initialized)
 *   rctx: region context (must be initialized)
 *
 * Returns: 0 on success, -1 on error
 */
int hotpants_setup_spatial_convolve(const config_struct* config,
                                    const kernel_context_struct* kctx,
                                    region_context_struct* rctx) {
  if (!config || !kctx || !rctx) return -1;

  set_kernel_context_globals(kctx, config);
  broadcast_region_context_globals(rctx);

  return 0;
}

/* =====================================================================
 * Context-Aware Core Algorithm Implementations (v2)
 * =====================================================================
 */

int fitKernel_v2(stamp_struct* stamps, int n_stamps,
                 float* imRef, float* imConv, float* imNoise,
                 const config_struct* config,
                 const kernel_context_struct* kctx,
                 region_context_struct* rctx,
                 double* kernel_coeffs,
                 double* meansig, double* scatter, int* n_skipped) {
  if (!stamps || !imRef || !imConv || !config || !kctx || !rctx ||
      !kernel_coeffs || !meansig || !scatter || !n_skipped) {
    return -1;
  }

  /* Set up global state from contexts */
  if (hotpants_setup_fitkernel(config, kctx, rctx) < 0) {
    return -1;
  }

  /* Call the existing fitKernel implementation */
  fitKernel(stamps, imRef, imConv, imNoise, kernel_coeffs, meansig, scatter,
            n_skipped);

  return 0;
}

int spatial_convolve_v2(float* image, float** var_image,
                        int ny, int nx,
                        double* kernel_coeffs,
                        const config_struct* config,
                        const kernel_context_struct* kctx,
                        region_context_struct* rctx,
                        float* output, int* conv_mask) {
  int conv_method = 1;  /* Use FFT by default (FFT is mandatory) */

  if (!image || !kernel_coeffs || !config || !kctx || !rctx || !output ||
      ny <= 0 || nx <= 0) {
    return -1;
  }

  /* Set up global state from contexts */
  if (hotpants_setup_spatial_convolve(config, kctx, rctx) < 0) {
    return -1;
  }

  /* Copy kernel coefficients into threadprivate global buffer */
  memcpy(rctx->kernel_coeffs, kernel_coeffs,
         (kctx->basis.n_comp_total + 1) * sizeof(double));

  /* Call the existing spatial_convolve implementation */
  spatial_convolve(image, var_image, ny, nx, rctx->kernel_coeffs, output,
                   &conv_method);

  /* If conv_mask was requested, populate it from mRData */
  if (conv_mask) {
    int pix_count = ny * nx;
    for (int i = 0; i < pix_count; i++) {
      conv_mask[i] = mRData ? mRData[i] : 0;
    }
  }

  return 0;
}
