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

/* Declare globals without TLS - we'll access them via proper extern with threadprivate pragma */
#ifdef _OPENMP
extern __thread int rPixX, rPixY;
extern __thread int* mRData;
#else
extern int rPixX, rPixY;
extern int* mRData;
#endif

extern char* forceConvolve;
extern int fwKernel, fwStamp;

/* Forward declarations - defined in functions.c */
extern void buildStamps(int sXMin, int sXMax, int sYMin, int sYMax, int* niS, int* ntS,
                        int getCenters, int rXBMin, int rYBMin, stamp_struct* ciStamps,
                        stamp_struct* ctStamps, float* iRData, float* tRData,
                        float hardX, float hardY);
extern void cutRegion(float* inArray, float* outArray, int xBMin, int xBMax,
                      int yBMin, int yBMax, int nx);
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
  int* region_i_mask_buffer;  /* Per-region science image mask */
  int* region_t_mask_buffer;  /* Per-region template mask */

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
    return -1; /* Already initialized */
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
  wrapper_context.region_i_mask_buffer = (int*)calloc(max_region_size, sizeof(int));
  wrapper_context.region_t_mask_buffer = (int*)calloc(max_region_size, sizeof(int));

  if (!wrapper_context.region_t_buffer || !wrapper_context.region_i_buffer ||
      !wrapper_context.region_mask_buffer || !wrapper_context.region_i_mask_buffer ||
      !wrapper_context.region_t_mask_buffer) {
    fprintf(stderr, "ERROR: Failed to allocate region buffers\n");
    return -1; /* Allocation failed */
  }

  /* Set global forceConvolve to "i" (convolve science image, matching template) */
  forceConvolve = "i";

  wrapper_context.initialized = 1;
  fprintf(stderr, "DEBUG: initBuildStampsContext succeeded\n");
  return 0;
}

/* Build stamps for a single region */
int buildStampsRegion(int region_x, int region_y,
                      stamp_struct* stamps, int n_stamps) {
  if (!wrapper_context.initialized) {
    return -1; /* Context not initialized */
  }

  /* Compute region boundaries */
  int rx_min, rx_max, ry_min, ry_max;
  int r_pix_x, r_pix_y;
  get_region_bounds(region_x, region_y,
                    wrapper_context.n_regions_x, wrapper_context.n_regions_y,
                    wrapper_context.full_nx, wrapper_context.full_ny,
                    &rx_min, &rx_max, &ry_min, &ry_max,
                    &r_pix_x, &r_pix_y);

  /* Set global region dimensions (used by buildStamps internally) */
  rPixX = r_pix_x;
  rPixY = r_pix_y;

  /* Extract region data from full images using numpy slicing equivalent */
  /* This is done in Python, so here we just use the pre-extracted buffers */
  for (int y = 0; y < r_pix_y; y++) {
    for (int x = 0; x < r_pix_x; x++) {
      int region_idx = y * r_pix_x + x;
      int full_idx = (ry_min + y) * wrapper_context.full_nx + (rx_min + x);
      wrapper_context.region_t_buffer[region_idx] = wrapper_context.template_image[full_idx];
      wrapper_context.region_i_buffer[region_idx] = wrapper_context.science_image[full_idx];
    }
  }

  /* Initialize mask arrays */
  memset(wrapper_context.region_mask_buffer, 0, r_pix_x * r_pix_y * sizeof(int));
  memset(wrapper_context.region_i_mask_buffer, 0, r_pix_x * r_pix_y * sizeof(int));
  memset(wrapper_context.region_t_mask_buffer, 0, r_pix_x * r_pix_y * sizeof(int));

  /* Set the global mask array pointer */
  mRData = wrapper_context.region_mask_buffer;

  /* Set up mask from image data */
  makeInputMask(wrapper_context.region_i_buffer, wrapper_context.region_t_buffer,
                wrapper_context.region_mask_buffer);

  /* Build stamps - iterate through stamp grid */
  int niS = 0, ntS = 0;  /* Stamp counters for this region */
  int stamps_built = 0;

  for (int stamp_y = 0; stamp_y < wrapper_context.stamps_per_region_y; stamp_y++) {
    for (int stamp_x = 0; stamp_x < wrapper_context.stamps_per_region_x; stamp_x++) {
      if (stamps_built >= n_stamps) {
        break; /* Don't exceed expected stamp count */
      }

      /* Compute stamp boundaries in full image coordinates */
      int sXMin = rx_min + (stamp_x * r_pix_x) / wrapper_context.stamps_per_region_x;
      int sXMax = rx_min + (((stamp_x + 1) * r_pix_x) / wrapper_context.stamps_per_region_x) - 1;
      int sYMin = ry_min + (stamp_y * r_pix_y) / wrapper_context.stamps_per_region_y;
      int sYMax = ry_min + (((stamp_y + 1) * r_pix_y) / wrapper_context.stamps_per_region_y) - 1;

      /* Clamp to region boundaries */
      if (sXMax >= wrapper_context.full_nx) sXMax = wrapper_context.full_nx - 1;
      if (sYMax >= wrapper_context.full_ny) sYMax = wrapper_context.full_ny - 1;

      /* Call buildStamps with region-relative coordinates */
      /* Note: buildStamps signature expects rXBMin, rYBMin for offset calculation */
      buildStamps(sXMin, sXMax, sYMin, sYMax, &niS, &ntS, 1,
                  rx_min, ry_min,
                  &stamps[stamps_built],  /* Single stamp per call */
                  &stamps[stamps_built],  /* Both point to same stamp for simplicity */
                  wrapper_context.region_i_buffer,
                  wrapper_context.region_t_buffer,
                  (float)fwKernel, (float)fwStamp);

      stamps_built++;
    }
  }

  return stamps_built; /* Return number of stamps built */
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
  if (wrapper_context.region_i_mask_buffer) {
    free(wrapper_context.region_i_mask_buffer);
    wrapper_context.region_i_mask_buffer = NULL;
  }
  if (wrapper_context.region_t_mask_buffer) {
    free(wrapper_context.region_t_mask_buffer);
    wrapper_context.region_t_mask_buffer = NULL;
  }

  /* Reset global pointers */
  mRData = NULL;

  /* Mark context as uninitialized */
  wrapper_context.initialized = 0;
}
