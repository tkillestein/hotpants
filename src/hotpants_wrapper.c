#include <string.h>
#include <stdlib.h>
#include <stdio.h>

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
extern char* forceConvolve;
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

  fprintf(stderr, "DEBUG: initKernelGlobals: nCompKer=%d nC=%d "
          "fwKSStamp=%d fwStamp=%d fwKernel=%d kcStep=%d\n",
          nCompKer, nC, fwKSStamp, fwStamp, fwKernel, kcStep);

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

  /* Estimate maximum region size for buffer allocation */
  int max_region_x = (nx + n_regions_x - 1) / n_regions_x + 1;
  int max_region_y = (ny + n_regions_y - 1) / n_regions_y + 1;
  int max_region_size = max_region_x * max_region_y;

  fprintf(stderr, "DEBUG: Allocating %d bytes for region buffers (max_region_size=%d)\n",
          max_region_size * (int)sizeof(float), max_region_size);

  /* Allocate region buffers */
  wrapper_context.region_t_buffer = (float*)calloc(max_region_size, sizeof(float));
  wrapper_context.region_i_buffer = (float*)calloc(max_region_size, sizeof(float));
  wrapper_context.region_mask_buffer = (int*)calloc(max_region_size, sizeof(int));

  if (!wrapper_context.region_t_buffer || !wrapper_context.region_i_buffer ||
      !wrapper_context.region_mask_buffer) {
    fprintf(stderr, "ERROR: Failed to allocate region buffers\n");
    return -1;
  }

  /* Set global forceConvolve to "i" (convolve science image, matching template) */
  forceConvolve = "i";

  wrapper_context.initialized = 1;
  fprintf(stderr, "DEBUG: initBuildStampsContext succeeded\n");
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

  /* Compute region boundaries */
  int rx_min, rx_max, ry_min, ry_max;
  int r_pix_x, r_pix_y;
  get_region_bounds(region_x, region_y,
                    wrapper_context.n_regions_x, wrapper_context.n_regions_y,
                    wrapper_context.full_nx, wrapper_context.full_ny,
                    &rx_min, &rx_max, &ry_min, &ry_max,
                    &r_pix_x, &r_pix_y);

  /* Set global region dimensions (required by C buildStamps) */
  rPixX = r_pix_x;
  rPixY = r_pix_y;

  fprintf(stderr, "DEBUG: buildStampsRegion(%d, %d): region_size=%dx%d\n",
          region_x, region_y, r_pix_x, r_pix_y);

  /* Extract region data from full images */
  for (int y = 0; y < r_pix_y; y++) {
    for (int x = 0; x < r_pix_x; x++) {
      int region_idx = y * r_pix_x + x;
      int full_idx = (ry_min + y) * wrapper_context.full_nx + (rx_min + x);
      wrapper_context.region_t_buffer[region_idx] = wrapper_context.template_image[full_idx];
      wrapper_context.region_i_buffer[region_idx] = wrapper_context.science_image[full_idx];
    }
  }

  /* Initialize mask array */
  memset(wrapper_context.region_mask_buffer, 0, r_pix_x * r_pix_y * sizeof(int));

  /* Set the global mask array pointer */
  mRData = wrapper_context.region_mask_buffer;

  /* Set up mask from image data */
  makeInputMask(wrapper_context.region_i_buffer, wrapper_context.region_t_buffer,
                wrapper_context.region_mask_buffer);

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
    fprintf(stderr, "ERROR: fillStampsForRegion: context not initialized\n");
    return -1;
  }
  int k;
  for (k = 0; k < n_stamps; k++) {
    fillStamp(&stamps[k],
              wrapper_context.region_i_buffer,
              wrapper_context.region_t_buffer);
  }
  return 0;
}

/* Clean up buildStamps context */
void cleanupBuildStampsContext(void) {
  if (!wrapper_context.initialized) {
    return;
  }

  fprintf(stderr, "DEBUG: cleanupBuildStampsContext called\n");

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
