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

/* Forward declarations - defined in functions.c */
extern void makeInputMask(float* iRData, float* tRData, int* mRData);

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
