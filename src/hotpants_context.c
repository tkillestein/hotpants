#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hotpants_context.h"
#include "defaults.h"
#include "allocate.h"

/* Forward declarations from alard.c */
extern void getKernelVec(void);

/* Computed component counts from kernel configuration */
static void compute_component_counts(const config_struct* config,
                                     gaussian_basis_struct* basis) {
  int i;

  /* Default 3-Gaussian basis */
  basis->ngauss = 3;
  basis->deg_fixe = (int*)calloc(basis->ngauss, sizeof(int));
  basis->sigma_gauss = (float*)calloc(basis->ngauss, sizeof(float));

  if (!basis->deg_fixe || !basis->sigma_gauss) {
    free(basis->deg_fixe);
    free(basis->sigma_gauss);
    return;
  }

  /* Default degrees and sigmas (matching vargs.c defaults) */
  basis->deg_fixe[0] = 6;
  basis->deg_fixe[1] = 4;
  basis->deg_fixe[2] = 2;

  /* sigma_gauss stored as 1/(2*sigma^2) convention */
  basis->sigma_gauss[0] = 1.0f / (2.0f * 0.7f * 0.7f);
  basis->sigma_gauss[1] = 1.0f / (2.0f * 1.5f * 1.5f);
  basis->sigma_gauss[2] = 1.0f / (2.0f * 3.0f * 3.0f);

  /* Compute total kernel components from Gaussian basis */
  basis->n_comp_ker = 0;
  for (i = 0; i < basis->ngauss; i++)
    basis->n_comp_ker += ((basis->deg_fixe[i] + 1) * (basis->deg_fixe[i] + 2)) / 2;

  /* Spatial polynomial terms for kernel */
  basis->n_comp = ((config->kernel_order + 1) * (config->kernel_order + 2)) / 2;
  basis->n_c = basis->n_comp_ker + 2;

  /* Background components */
  basis->n_comp_bg = (basis->n_comp_ker - 1) * basis->n_comp + 1;
  basis->n_bg_vectors = ((config->bg_order + 1) * (config->bg_order + 2)) / 2;
  basis->n_comp_total = basis->n_comp_ker * basis->n_comp + basis->n_bg_vectors;
}

/* Helper: compute fwStamp from image and region dimensions */
static int compute_fw_stamp(int image_nx, int image_ny,
                            int n_regions_x, int n_regions_y,
                            int stamps_per_region_x, int stamps_per_region_y,
                            int fw_kernel, int fw_ks_stamp) {
  int min_dim_per_region_x = image_nx / n_regions_x / stamps_per_region_x;
  int min_dim_per_region_y = image_ny / n_regions_y / stamps_per_region_y;
  int min_dim = (min_dim_per_region_x < min_dim_per_region_y)
                ? min_dim_per_region_x : min_dim_per_region_y;

  int fw_stamp = min_dim - fw_kernel;
  if (fw_stamp % 2 == 0) fw_stamp -= 1;  /* keep odd */

  /* Sanity check */
  if (fw_stamp < fw_ks_stamp) {
    fw_stamp = fw_ks_stamp + fw_kernel;
    if (fw_stamp % 2 == 0) fw_stamp -= 1;
  }

  return fw_stamp;
}

/* Helper: allocate kernel basis vectors via getKernelVec() callback */
static int allocate_kernel_basis(kernel_context_struct* ctx) {
  int i;
  int fw_kernel = ctx->fw_kernel;
  int n_comp_ker = ctx->basis.n_comp_ker;

  /* Allocate kernel_vec[][] */
  ctx->basis.kernel_vec = (double**)calloc(n_comp_ker, sizeof(double*));
  if (!ctx->basis.kernel_vec) return -1;

  for (i = 0; i < n_comp_ker; i++) {
    ctx->basis.kernel_vec[i] = (double*)calloc(fw_kernel * fw_kernel, sizeof(double));
    if (!ctx->basis.kernel_vec[i]) {
      while (i >= 0) {
        free(ctx->basis.kernel_vec[i]);
        i--;
      }
      free(ctx->basis.kernel_vec);
      ctx->basis.kernel_vec = NULL;
      return -1;
    }
  }

  /* Allocate filter_x and filter_y */
  ctx->basis.filter_x = (double*)calloc(n_comp_ker * fw_kernel, sizeof(double));
  ctx->basis.filter_y = (double*)calloc(n_comp_ker * fw_kernel, sizeof(double));
  if (!ctx->basis.filter_x || !ctx->basis.filter_y) {
    free(ctx->basis.filter_x);
    free(ctx->basis.filter_y);
    return -1;
  }

  /* Call getKernelVec() to fill the basis (requires some globals set) */
  /* NOTE: getKernelVec() reads deg_fixe, sigma_gauss, ngauss globals.
   * We have set these already; getKernelVec() will populate kernel_vec, filter_x/y. */
  getKernelVec();

  return 0;
}

int hotpants_create_kernel_context(const config_struct* config,
                                   int image_nx, int image_ny,
                                   int n_regions_x, int n_regions_y,
                                   int stamps_per_region_x, int stamps_per_region_y,
                                   kernel_context_struct** out_ctx) {
  kernel_context_struct* ctx = NULL;

  if (!config || !out_ctx) return -1;

  ctx = (kernel_context_struct*)calloc(1, sizeof(kernel_context_struct));
  if (!ctx) return -1;

  /* Compute frame widths */
  ctx->fw_kernel = 2 * config->kernel_half_width + 1;
  ctx->fw_ks_stamp = 2 * config->hw_ks_stamp + 1;
  ctx->kc_step = ctx->fw_kernel;

  /* Compute fwStamp from image and region layout */
  ctx->fw_stamp = compute_fw_stamp(image_nx, image_ny,
                                   n_regions_x, n_regions_y,
                                   stamps_per_region_x, stamps_per_region_y,
                                   ctx->fw_kernel, ctx->fw_ks_stamp);

  /* Compute component counts and allocate Gaussian basis arrays */
  compute_component_counts(config, &ctx->basis);
  if (!ctx->basis.deg_fixe || !ctx->basis.sigma_gauss) {
    free(ctx);
    return -1;
  }

  /* Allocate kernel basis vectors (kernel_vec, filter_x/y) */
  if (allocate_kernel_basis(ctx) < 0) {
    free(ctx->basis.deg_fixe);
    free(ctx->basis.sigma_gauss);
    free(ctx);
    return -1;
  }

  ctx->is_initialized = 1;
  *out_ctx = ctx;
  return 0;
}

void hotpants_cleanup_kernel_context(kernel_context_struct* ctx) {
  int i;
  if (!ctx) return;

  if (ctx->basis.kernel_vec) {
    for (i = 0; i < ctx->basis.n_comp_ker; i++) {
      if (ctx->basis.kernel_vec[i]) free(ctx->basis.kernel_vec[i]);
    }
    free(ctx->basis.kernel_vec);
    ctx->basis.kernel_vec = NULL;
  }

  if (ctx->basis.filter_x) {
    free(ctx->basis.filter_x);
    ctx->basis.filter_x = NULL;
  }

  if (ctx->basis.filter_y) {
    free(ctx->basis.filter_y);
    ctx->basis.filter_y = NULL;
  }

  if (ctx->basis.deg_fixe) {
    free(ctx->basis.deg_fixe);
    ctx->basis.deg_fixe = NULL;
  }

  if (ctx->basis.sigma_gauss) {
    free(ctx->basis.sigma_gauss);
    ctx->basis.sigma_gauss = NULL;
  }

  free(ctx);
}

int hotpants_create_region_context(const kernel_context_struct* kctx,
                                   const config_struct* config,
                                   int n_stamps,
                                   region_context_struct** out_ctx) {
  region_context_struct* ctx = NULL;
  int i;

  if (!kctx || !config || !out_ctx || n_stamps <= 0) return -1;

  ctx = (region_context_struct*)calloc(1, sizeof(region_context_struct));
  if (!ctx) return -1;

  /* Allocate temporary buffers */
  ctx->temp = (float*)calloc((kctx->fw_ks_stamp + kctx->fw_kernel) * kctx->fw_ks_stamp,
                             sizeof(float));
  ctx->temp2 = (float*)calloc((kctx->fw_ks_stamp + kctx->fw_kernel) * kctx->fw_ks_stamp,
                              sizeof(float));

  /* Allocate algorithm working buffers */
  ctx->kernel = (double*)calloc(kctx->fw_kernel * kctx->fw_kernel, sizeof(double));
  ctx->kernel_coeffs = (double*)calloc(kctx->basis.n_comp_total + 1, sizeof(double));
  ctx->check_stack = (double*)calloc(n_stamps * kctx->basis.n_comp, sizeof(double));

  /* Weight and check matrices */
  ctx->wxy = (double**)calloc(n_stamps, sizeof(double*));
  ctx->check_mat = (double**)calloc(n_stamps, sizeof(double*));

  if (!ctx->temp || !ctx->temp2 || !ctx->kernel || !ctx->kernel_coeffs ||
      !ctx->check_stack || !ctx->wxy || !ctx->check_mat) {
    hotpants_cleanup_region_context(ctx);
    return -1;
  }

  /* Allocate rows for wxy and check_mat */
  for (i = 0; i < n_stamps; i++) {
    ctx->wxy[i] = (double*)calloc(kctx->basis.n_comp, sizeof(double));
    ctx->check_mat[i] = (double*)calloc(kctx->basis.n_comp, sizeof(double));

    if (!ctx->wxy[i] || !ctx->check_mat[i]) {
      hotpants_cleanup_region_context(ctx);
      return -1;
    }
  }

  /* Allocate check_vec */
  ctx->check_vec = (double*)calloc(kctx->basis.n_comp, sizeof(double));
  if (!ctx->check_vec) {
    hotpants_cleanup_region_context(ctx);
    return -1;
  }

  ctx->n_stamps_allocated = n_stamps;
  ctx->is_initialized = 1;
  *out_ctx = ctx;
  return 0;
}

void hotpants_cleanup_region_context(region_context_struct* ctx) {
  int i;
  if (!ctx) return;

  if (ctx->wxy) {
    for (i = 0; i < ctx->n_stamps_allocated; i++) {
      if (ctx->wxy[i]) free(ctx->wxy[i]);
    }
    free(ctx->wxy);
    ctx->wxy = NULL;
  }

  if (ctx->check_mat) {
    for (i = 0; i < ctx->n_stamps_allocated; i++) {
      if (ctx->check_mat[i]) free(ctx->check_mat[i]);
    }
    free(ctx->check_mat);
    ctx->check_mat = NULL;
  }

  if (ctx->temp) {
    free(ctx->temp);
    ctx->temp = NULL;
  }

  if (ctx->temp2) {
    free(ctx->temp2);
    ctx->temp2 = NULL;
  }

  if (ctx->kernel) {
    free(ctx->kernel);
    ctx->kernel = NULL;
  }

  if (ctx->kernel_coeffs) {
    free(ctx->kernel_coeffs);
    ctx->kernel_coeffs = NULL;
  }

  if (ctx->check_stack) {
    free(ctx->check_stack);
    ctx->check_stack = NULL;
  }

  if (ctx->check_vec) {
    free(ctx->check_vec);
    ctx->check_vec = NULL;
  }

  free(ctx);
}

int hotpants_default_config(config_struct* out_config) {
  if (!out_config) return -1;

  memset(out_config, 0, sizeof(config_struct));

  /* Kernel configuration defaults */
  out_config->kernel_half_width = D_HWKERNEL;
  out_config->kernel_order = D_KORDER;
  out_config->bg_order = D_BGORDER;
  out_config->fit_threshold = D_KFITTHRESH;
  out_config->scale_fit_threshold = D_SFITTHRESH;
  out_config->n_ks_stamps = D_NKSSTAMPS;
  out_config->hw_ks_stamp = D_HWKSSTAMP;

  /* Thresholds (neutral: no clipping) */
  out_config->template_upper_threshold = D_UTHRESH;
  out_config->template_lower_threshold = D_LTHRESH;
  out_config->science_upper_threshold = D_UTHRESH;
  out_config->science_lower_threshold = D_LTHRESH;

  /* Gains and noise (unity) */
  out_config->template_gain = D_GAIN;
  out_config->science_gain = D_GAIN;
  out_config->template_readnoise = D_RDNOISE;
  out_config->science_readnoise = D_RDNOISE;
  out_config->template_pedestal = D_PEDESTAL;
  out_config->science_pedestal = D_PEDESTAL;

  /* Tuning parameters */
  out_config->stat_sig = D_STATSIG;
  out_config->ker_sig_reject = D_KSIGREJECT;
  out_config->ker_frac_mask = D_KFRACMASK;
  out_config->fill_val = D_FILL;
  out_config->fill_val_noise = D_FILLNOISE;
  out_config->kf_spread_mask1 = D_INMASKFSPREAD;
  out_config->kf_spread_mask2 = D_OUMASKFSPREAD;
  out_config->min_frac_good_stamps = D_NFITTHRESH;

  /* Region/stamp layout */
  out_config->n_regions_x = 1;
  out_config->n_regions_y = 1;
  out_config->stamps_per_region_x = D_NSTAMPS;
  out_config->stamps_per_region_y = D_NSTAMPS;
  out_config->use_full_ss = D_USEFULLSS;

  /* Runtime options */
  out_config->verbose = D_VERBOSE;
  out_config->n_thread = D_NTHREAD;

  /* Direction: template ⊗ K ≈ science (standard) */
  out_config->force_convolve = 't';

  return 0;
}
