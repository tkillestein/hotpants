#ifndef HOTPANTS_CONTEXT_H
#define HOTPANTS_CONTEXT_H

/*
 * HOTPANTS Configuration & Context Structures
 *
 * This header defines data structures that encapsulate global state,
 * eliminating the need for global variables. Structures are passed explicitly
 * through the call chain, improving thread safety and code clarity.
 *
 * Three main structures:
 * - config_struct: User-settable configuration (kernel, noise, fit parameters)
 * - kernel_context_struct: Pre-computed derived values & Gaussian basis
 * - region_context_struct: Per-region working buffers (thread-local)
 *
 * Reference: Alard & Lupton (1998), ApJ 503:325
 */

#include <stddef.h>

/* =====================================================================
 * Configuration Structure
 * =====================================================================
 * Holds all user-settable parameters for kernel fitting and convolution.
 * Passed to initialization functions and wrapped in context structures.
 */
typedef struct {
  /* Kernel basis configuration */
  int kernel_half_width;           /* hwKernel: half-width of kernel region (pixels) */
  int kernel_order;                /* kerOrder: spatial polynomial order */
  int bg_order;                    /* bgOrder: background polynomial order */
  int n_ks_stamps;                 /* nKSStamps: kernel substamps per stamp */
  int hw_ks_stamp;                 /* hwKSStamp: half-width of kernel substamp */

  /* Noise/threshold configuration */
  float template_upper_threshold;  /* tUThresh: template saturation threshold */
  float template_lower_threshold;  /* tLThresh: template lower valid threshold */
  float science_upper_threshold;   /* iUThresh: science saturation threshold */
  float science_lower_threshold;   /* iLThresh: science lower valid threshold */
  float template_gain;             /* tGain: template gain (e-/ADU) */
  float science_gain;              /* iGain: science gain (e-/ADU) */
  float template_readnoise;        /* tRdnoise: template readnoise (e-) */
  float science_readnoise;         /* iRdnoise: science readnoise (e-) */
  float template_pedestal;         /* tPedestal: template pedestal (ADU) */
  float science_pedestal;          /* iPedestal: science pedestal (ADU) */

  /* Fit configuration */
  float fit_threshold;             /* kerFitThresh: sigma threshold for fit */
  float scale_fit_threshold;       /* scaleFitThresh: scale factor for retry threshold */
  float min_frac_good_stamps;      /* minFracGoodStamps: minimum fraction of valid stamps */

  /* Tuning parameters */
  float stat_sig;                  /* statSig: sigma-clipping threshold for statistics */
  float ker_sig_reject;            /* kerSigReject: sigma rejection for stamp check */
  float ker_frac_mask;             /* kerFracMask: mask fraction threshold */
  float fill_val;                  /* fillVal: fill value for bad pixels */
  float fill_val_noise;            /* fillValNoise: fill value for noise image */
  float kf_spread_mask1;           /* kfSpreadMask1: kernel fit spread mask threshold 1 */
  float kf_spread_mask2;           /* kfSpreadMask2: kernel fit spread mask threshold 2 */

  /* Region/stamp layout */
  int n_regions_x;                 /* nRegX: regions per axis X */
  int n_regions_y;                 /* nRegY: regions per axis Y */
  int stamps_per_region_x;         /* nStampX: stamps per region X */
  int stamps_per_region_y;         /* nStampY: stamps per region Y */
  int use_full_ss;                 /* useFullSS: use full substamp division */

  /* Runtime options */
  int verbose;                     /* verbose: verbosity level */
  int n_thread;                    /* nThread: number of OpenMP threads */

  /* Direction: 't' = template ⊗ K ≈ science, 's' = science ⊗ K ≈ template */
  char force_convolve;
} config_struct;

/* =====================================================================
 * Gaussian Basis Structure (part of kernel context)
 * =====================================================================
 * Holds precomputed Gaussian basis parameters and basis vectors.
 */
typedef struct {
  int ngauss;                      /* Number of Gaussian basis components */
  int* deg_fixe;                   /* [ngauss] Polynomial degree for each Gaussian */
  float* sigma_gauss;              /* [ngauss] Sigma parameters (1/(2*sigma^2) convention) */

  int n_comp_ker;                  /* nCompKer: total components in Gaussian basis */
  int n_comp;                      /* nComp: polynomial terms for kernel spatial variation */
  int n_c;                         /* nC: used by allocateStamps (nCompKer + 2) */
  int n_comp_bg;                   /* nCompBG: background components */
  int n_bg_vectors;                /* nBGVectors: background polynomial terms */
  int n_comp_total;                /* nCompTotal: total components (kernel + background) */

  /* Kernel basis vectors (precomputed once per config) */
  double** kernel_vec;             /* [n_comp_ker][fw_kernel*fw_kernel] */
  double* filter_x;                /* [n_comp_ker * fw_kernel] separable x-filter */
  double* filter_y;                /* [n_comp_ker * fw_kernel] separable y-filter */
} gaussian_basis_struct;

/* =====================================================================
 * Kernel Context Structure
 * =====================================================================
 * Holds all derived values computed from config_struct, including
 * frame widths, component counts, and precomputed Gaussian basis.
 * Allocated once per configuration; valid across multiple regions.
 */
typedef struct {
  /* Derived configuration (computed from config_struct) */
  int fw_kernel;                   /* fwKernel = 2*kernel_half_width + 1 */
  int fw_stamp;                    /* fwStamp: stamp frame width (region-dependent) */
  int fw_ks_stamp;                 /* fwKSStamp = 2*hw_ks_stamp + 1 */
  int kc_step;                     /* kcStep: convolution step size (= fw_kernel) */

  /* Gaussian basis and derived counts */
  gaussian_basis_struct basis;

  /* Session state for stamp filling */
  int is_initialized;              /* Flag: context valid and allocated */
} kernel_context_struct;

/* =====================================================================
 * Region Context Structure
 * =====================================================================
 * Holds all per-region working buffers and state.
 * Allocated per region; thread-local when using OpenMP.
 * Mirrors the threadprivate globals from globals.h.
 */
typedef struct {
  /* Per-region dimensions */
  int pix_x;                       /* rPixX: region width in pixels */
  int pix_y;                       /* rPixY: region height in pixels */
  int* mask;                       /* mRData: region mask array [pix_y*pix_x] */

  /* Scratch buffers for stamp processing */
  float* temp;                     /* Temporary buffer for make_model/getStampSig */
  float* temp2;                    /* Additional temporary buffer */

  /* Algorithm working buffers (dimensions depend on kernel context) */
  double** wxy;                    /* Weight matrix [n_stamps][n_comp] */
  double* kernel;                  /* Kernel evaluation buffer [fw_kernel*fw_kernel] */
  double* kernel_coeffs;           /* Kernel coefficients [n_comp_total+1] */
  double* check_stack;             /* Sigma-clipping stack */
  double** check_mat;              /* Check matrix [n_stamps][n_comp] */
  double* check_vec;               /* Check vector [n_comp] */

  /* Allocation tracking */
  int n_stamps_allocated;          /* Number of stamps allocated for this region */
  int is_initialized;              /* Flag: context valid and allocated */
} region_context_struct;

/* =====================================================================
 * Public API: Context Initialization & Cleanup
 * =====================================================================
 */

/*
 * Create and initialize a kernel context from configuration.
 *
 * Allocates derived values (frame widths, component counts) and precomputes
 * the Gaussian basis vectors (kernel_vec, filter_x/y). This context is
 * valid across multiple regions and should be reused for efficiency.
 *
 * Args:
 *   config: configuration struct (read-only)
 *   image_nx, image_ny: full image dimensions (used to compute fwStamp)
 *   n_regions_x, n_regions_y: region grid (used to compute fwStamp)
 *   stamps_per_region_x, stamps_per_region_y: stamp grid (used to compute fwStamp)
 *   out_ctx: pointer to kernel_context_struct* (filled on success)
 *
 * Returns:
 *   0 on success, -1 on allocation failure
 *
 * Allocates: kernel_vec[], filter_x[], filter_y[], deg_fixe[], sigma_gauss[]
 * Must be freed with hotpants_cleanup_kernel_context().
 */
int hotpants_create_kernel_context(const config_struct* config,
                                   int image_nx, int image_ny,
                                   int n_regions_x, int n_regions_y,
                                   int stamps_per_region_x, int stamps_per_region_y,
                                   kernel_context_struct** out_ctx);

/*
 * Free all memory associated with a kernel context.
 *
 * Args:
 *   ctx: kernel context to free (may be NULL)
 */
void hotpants_cleanup_kernel_context(kernel_context_struct* ctx);

/*
 * Create and initialize a region context from kernel context.
 *
 * Allocates per-region working buffers: temp, temp2, wxy, kernel, kernel_coeffs,
 * check_stack, check_mat, check_vec. Buffer sizes depend on kernel context.
 *
 * Args:
 *   kctx: kernel context (read-only; must be valid)
 *   config: configuration struct (read-only; used for derived size calculations)
 *   n_stamps: number of stamps for this region
 *   out_ctx: pointer to region_context_struct* (filled on success)
 *
 * Returns:
 *   0 on success, -1 on allocation failure
 *
 * Allocates: temp, temp2, wxy, kernel, kernel_coeffs, check_stack, check_mat, check_vec
 * Must be freed with hotpants_cleanup_region_context().
 */
int hotpants_create_region_context(const kernel_context_struct* kctx,
                                   const config_struct* config,
                                   int n_stamps,
                                   region_context_struct** out_ctx);

/*
 * Free all memory associated with a region context.
 *
 * Args:
 *   ctx: region context to free (may be NULL)
 */
void hotpants_cleanup_region_context(region_context_struct* ctx);

/*
 * Create a default configuration struct with standard parameters.
 *
 * Fills config_struct with sensible defaults matching the HOTPANTS CLI defaults:
 *   kernel_half_width = 15
 *   kernel_order = 2
 *   bg_order = 1
 *   fit_threshold = 20.0
 *   scale_fit_threshold = 0.5
 *   n_ks_stamps = 3
 *   hw_ks_stamp = 10
 *   (thresholds and gains set to neutral/unity values)
 *
 * Args:
 *   out_config: pointer to config_struct (filled on return)
 *
 * Returns:
 *   0 on success
 */
int hotpants_default_config(config_struct* out_config);

#endif /* HOTPANTS_CONTEXT_H */
