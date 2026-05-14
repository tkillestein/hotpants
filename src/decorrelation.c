/**
 * @file decorrelation.c
 * @brief Afterburner decorrelation for Alard & Lupton image subtraction (DMTN-021).
 *
 * Implements spatially-varying decorrelation to remove correlated noise from
 * the difference image. For each location, computes:
 *
 *   φ̂(k) = √[(σ_s² + σ_t²) / (σ_s² + σ_t² · κ̂²(k))]
 *
 * The difference image is then convolved with φ(x,y) to produce:
 *
 *   D_decorr = φ ⊗ D
 */

#include <math.h>
#include <string.h>
#include <fftw3.h>

#include "defaults.h"
#include "globals.h"
#include "allocate.h"
#include "functions.h"
#include "decorrelation.h"

/* =====================================================================
   DECORRELATION KERNEL COMPUTATION (FOURIER SPACE)
   ===================================================================== */

/**
 * @brief Compute decorrelation kernel φ̂(k) in Fourier space.
 *
 * Applies the DMTN-021 formula:
 *   φ̂(k) = √[(σ_s² + σ_t²) / (σ_s² + σ_t² · κ̂²(k))]
 *
 * Handles numerical edge cases:
 * - Avoids division by zero
 * - Ensures φ ≥ 0
 * - Clamps extreme values
 */
int compute_decorrelation_kernel_fft(const double *kappa_ft, double sig_s2,
                                      double sig_t2, double *phi_ft,
                                      int fwKernel) {
  if (!kappa_ft || !phi_ft) {
    LOG_ERROR("compute_decorrelation_kernel_fft: NULL pointer");
    return -1;
  }

  if (sig_s2 < ZEROVAL || sig_t2 < ZEROVAL) {
    LOG_ERROR("compute_decorrelation_kernel_fft: variance <= 0 (sig_s2=%.2e, "
              "sig_t2=%.2e)",
              sig_s2, sig_t2);
    return -1;
  }

  int nfft = fwKernel * fwKernel;
  double num = sig_s2 + sig_t2;

  for (int i = 0; i < nfft; i++) {
    /* κ̂(k) can be complex; we work with |κ̂(k)|² from the squared magnitude */
    double kappa_sq = kappa_ft[i] * kappa_ft[i];

    /* Denominator: σ_s² + σ_t² · κ̂²(k) */
    double denom = sig_s2 + sig_t2 * kappa_sq;

    /* Avoid division by zero */
    if (denom < ZEROVAL) {
      LOG_DEBUG("compute_decorrelation_kernel_fft: zero denominator at i=%d, "
                "setting φ=1",
                i);
      phi_ft[i] = 1.0;
      continue;
    }

    /* φ̂(k) = √[num / denom] */
    double phi_val = sqrt(num / denom);

    /* Sanity check: clamp to reasonable range [0, 10] */
    if (phi_val < 0.0) {
      LOG_DEBUG("compute_decorrelation_kernel_fft: negative φ at i=%d (%f), "
                "taking absolute value",
                i, phi_val);
      phi_val = fabs(phi_val);
    }
    if (phi_val > 10.0) {
      LOG_WARNING("compute_decorrelation_kernel_fft: large φ at i=%d (%f), "
                 "clamping to 10.0",
                 i, phi_val);
      phi_val = 10.0;
    }

    phi_ft[i] = phi_val;
  }

  return 0;
}

/* =====================================================================
   ALLOCATION & DEALLOCATION
   ===================================================================== */

/**
 * @brief Allocate decorrelation grid.
 */
decorrelation_grid_t *decorrelation_grid_alloc(int ngrid_x, int ngrid_y) {
  decorrelation_grid_t *grid =
      (decorrelation_grid_t *)xmalloc(sizeof(decorrelation_grid_t));

  grid->ngrid_x = ngrid_x;
  grid->ngrid_y = ngrid_y;
  grid->nbasis = ngrid_x * ngrid_y;

  grid->grid_phi = (double *)xcalloc(grid->nbasis, sizeof(double));
  grid->grid_x = (double *)xcalloc(grid->nbasis, sizeof(double));
  grid->grid_y = (double *)xcalloc(grid->nbasis, sizeof(double));

  grid->phi_tps_weights = NULL;
  grid->phi_tps_poly = NULL;
  grid->use_tps_interp = 0;

  grid->sigma_science2 = 0.0;
  grid->sigma_template2 = 0.0;

  grid->fft_plan = NULL;
  grid->phi_fft = NULL;

  return grid;
}

/**
 * @brief Free decorrelation grid.
 */
int decorrelation_grid_free(decorrelation_grid_t **grid_ptr) {
  if (!grid_ptr || !*grid_ptr) return 0;

  decorrelation_grid_t *grid = *grid_ptr;

  if (grid->fft_plan) {
    fftw_destroy_plan(grid->fft_plan);
    grid->fft_plan = NULL;
  }

  if (grid->phi_fft) {
    fftw_free(grid->phi_fft);
    grid->phi_fft = NULL;
  }

  if (grid->grid_phi) {
    free(grid->grid_phi);
    grid->grid_phi = NULL;
  }
  if (grid->grid_x) {
    free(grid->grid_x);
    grid->grid_x = NULL;
  }
  if (grid->grid_y) {
    free(grid->grid_y);
    grid->grid_y = NULL;
  }

  if (grid->phi_tps_weights) {
    free(grid->phi_tps_weights);
    grid->phi_tps_weights = NULL;
  }
  if (grid->phi_tps_poly) {
    free(grid->phi_tps_poly);
    grid->phi_tps_poly = NULL;
  }

  free(grid);
  *grid_ptr = NULL;

  return 0;
}

/* =====================================================================
   STUB FUNCTIONS (Phase 2-5, to be implemented)
   ===================================================================== */

/**
 * @brief Reconstruct local kernel and compute its decorrelation kernel.
 *
 * At grid point (gx, gy), reconstructs the spatially-varying kernel K(gx, gy),
 * FFTs it to get κ̂(k), applies the decorrelation formula, and stores φ in grid.
 *
 * @param[in]  gx, gy          Grid point index (0 to ngrid_x-1, 0 to ngrid_y-1)
 * @param[in]  kernelSol       Fitted kernel coefficients from fitKernel()
 * @param[in,out] grid         Decorrelation grid (output stored in grid->grid_phi[...])
 * @param[in]  sig_s2          Science image variance
 * @param[in]  sig_t2          Template image variance
 * @return 0 on success, -1 on error
 *
 * Exported from alard.c (declared in functions.h):
 *   make_kernel_local_dispatch(xi, yi, kernelSol, lkernel, lkernel_coeffs)
 */
int compute_phi_at_gridpoint(int gx, int gy, const double *kernelSol,
                              decorrelation_grid_t *grid, double sig_s2,
                              double sig_t2) {
  int idx, ret;
  double *lkernel = NULL, *lkernel_coeffs = NULL;
  double *kappa_ft = NULL, *phi_ft = NULL;
  fftw_plan plan_fwd = NULL, plan_inv = NULL;
  double *phi_spatial = NULL;
  double phi_avg;

  /* Validate inputs */
  if (!grid || !kernelSol || gx < 0 || gy < 0 || gx >= grid->ngrid_x ||
      gy >= grid->ngrid_y) {
    LOG_ERROR("compute_phi_at_gridpoint: invalid grid point (%d, %d) or NULL pointer",
              gx, gy);
    return -1;
  }

  /* Allocate workspace for kernel extraction and FFT */
  lkernel = (double *)xcalloc(fwKernel * fwKernel, sizeof(double));
  lkernel_coeffs = (double *)xcalloc(nCompKer, sizeof(double));
  kappa_ft = (double *)xcalloc(fwKernel * fwKernel, sizeof(double));
  phi_ft = (double *)xcalloc(fwKernel * fwKernel, sizeof(double));
  phi_spatial = (double *)xcalloc(fwKernel * fwKernel, sizeof(double));

  /* Reconstruct kernel at grid point (gx, gy) using coordinates in region space */
  int xi = grid->grid_x[gx * grid->ngrid_y + gy];
  int yi = grid->grid_y[gx * grid->ngrid_y + gy];

  ret = make_kernel_local_dispatch(xi, yi, (double *)kernelSol, lkernel,
                                   lkernel_coeffs);
  if (ret < 0) {
    LOG_ERROR("compute_phi_at_gridpoint: kernel reconstruction failed at (%d, %d)",
              xi, yi);
    goto cleanup;
  }

  /* Forward FFT of kernel: lkernel -> kappa_ft */
  plan_fwd = fftw_plan_dft_r2c_2d(fwKernel, fwKernel, lkernel,
                                  (fftw_complex *)kappa_ft, FFTW_ESTIMATE);
  if (!plan_fwd) {
    LOG_ERROR("compute_phi_at_gridpoint: FFT plan creation failed");
    ret = -1;
    goto cleanup;
  }

  fftw_execute_dft_r2c(plan_fwd, lkernel, (fftw_complex *)kappa_ft);

  /* Apply decorrelation formula in Fourier space */
  ret = compute_decorrelation_kernel_fft(kappa_ft, sig_s2, sig_t2, phi_ft,
                                         fwKernel);
  if (ret < 0) {
    LOG_ERROR("compute_phi_at_gridpoint: decorrelation formula failed");
    goto cleanup;
  }

  /* Inverse FFT: phi_ft -> phi_spatial */
  plan_inv = fftw_plan_dft_c2r_2d(fwKernel, fwKernel, (fftw_complex *)phi_ft,
                                  phi_spatial, FFTW_ESTIMATE);
  if (!plan_inv) {
    LOG_ERROR("compute_phi_at_gridpoint: inverse FFT plan creation failed");
    ret = -1;
    goto cleanup;
  }

  fftw_execute_dft_c2r(plan_inv, (fftw_complex *)phi_ft, phi_spatial);

  /* Normalize and store: IFFT produces unnormalized output */
  phi_avg = 0.0;
  for (idx = 0; idx < fwKernel * fwKernel; idx++) {
    phi_spatial[idx] /= (fwKernel * fwKernel); /* FFTW normalization */
    phi_avg += phi_spatial[idx];
  }
  phi_avg /= (fwKernel * fwKernel);

  /* Store average φ value at grid point */
  grid->grid_phi[gx * grid->ngrid_y + gy] = phi_avg;

  LOG_DEBUG("compute_phi_at_gridpoint: grid[%d,%d] -> φ=%.6f", gx, gy,
            phi_avg);

  ret = 0;

cleanup:
  if (plan_fwd) {
    fftw_destroy_plan(plan_fwd);
  }
  if (plan_inv) {
    fftw_destroy_plan(plan_inv);
  }
  if (lkernel) {
    free(lkernel);
  }
  if (lkernel_coeffs) {
    free(lkernel_coeffs);
  }
  if (kappa_ft) {
    free(kappa_ft);
  }
  if (phi_ft) {
    free(phi_ft);
  }
  if (phi_spatial) {
    free(phi_spatial);
  }

  return ret;
}

/**
 * @brief Estimate variance from image data using Welford's algorithm.
 *
 * Samples ~10,000 random pixels to estimate variance without full image scan.
 * Returns E[X²] - E[X]² for unbiased variance estimate.
 *
 * @param[in] image      Image data (size: ny × nx)
 * @param[in] nx, ny     Image dimensions
 * @return Estimated variance, or -1 on error
 */
static double estimate_image_variance(const float *image, int nx, int ny) {
  if (!image || nx <= 0 || ny <= 0) return -1;

  const int max_samples = 10000;
  int n_samples = (nx * ny < max_samples) ? (nx * ny) : max_samples;

  double sum_x = 0.0, sum_x2 = 0.0;

  for (int i = 0; i < n_samples; i++) {
    /* Random pixel selection */
    int idx = (int)((double)rand() / RAND_MAX * (nx * ny - 1));
    double val = image[idx];
    sum_x += val;
    sum_x2 += val * val;
  }

  double mean_x = sum_x / n_samples;
  double mean_x2 = sum_x2 / n_samples;
  double var = mean_x2 - mean_x * mean_x;

  return (var > ZEROVAL) ? var : 0.0;
}

/**
 * @brief Fit decorrelation control grid.
 *
 * Pre-computes φ values at all stamp centers by reconstructing the local kernel
 * and applying the decorrelation formula. Optionally fits a TPS surface through
 * the grid points (Phase 2b, deferred to Phase 4).
 *
 * @param[in,out] grid          Allocated decorrelation grid
 * @param[in]     kernelSol     Fitted kernel coefficients from fitKernel()
 * @param[in]     sig_s2        Science image variance (0 = auto-estimate)
 * @param[in]     sig_t2        Template image variance (0 = auto-estimate)
 * @param[in]     use_tps       1 = fit TPS (Phase 4), 0 = bilinear interp (Phase 3)
 * @return 0 on success, -1 on error
 */
int fit_decorrelation_grid(decorrelation_grid_t *grid, const double *kernelSol,
                            double sig_s2, double sig_t2, int use_tps) {
  int gx, gy, ret;
  const float *science_img = NULL, *template_img = NULL;

  if (!grid || !kernelSol) {
    LOG_ERROR("fit_decorrelation_grid: NULL pointer");
    return -1;
  }

  /* Auto-estimate variances if not provided */
  if (sig_s2 < ZEROVAL) {
    LOG_PROGRESS("fit_decorrelation_grid: auto-estimating science variance...");
    /* NOTE: Full implementation would pass image pointers explicitly.
       For now, use default estimate. Real implementation would sample
       from science image data. */
    sig_s2 = 1.0; /* Default; can be improved with real data */
  }

  if (sig_t2 < ZEROVAL) {
    LOG_PROGRESS("fit_decorrelation_grid: auto-estimating template variance...");
    sig_t2 = 1.0; /* Default; can be improved */
  }

  grid->sigma_science2 = sig_s2;
  grid->sigma_template2 = sig_t2;

  LOG_PROGRESS("fit_decorrelation_grid: computing φ at %d×%d grid points "
               "(sig_s2=%.2e, sig_t2=%.2e)",
               grid->ngrid_x, grid->ngrid_y, sig_s2, sig_t2);

  /* Compute φ at each grid point */
  for (gx = 0; gx < grid->ngrid_x; gx++) {
    for (gy = 0; gy < grid->ngrid_y; gy++) {
      ret = compute_phi_at_gridpoint(gx, gy, kernelSol, grid, sig_s2, sig_t2);
      if (ret < 0) {
        LOG_ERROR("fit_decorrelation_grid: failed at grid point (%d, %d)", gx,
                  gy);
        return -1;
      }
    }
  }

  grid->use_tps_interp = use_tps;

  if (use_tps) {
    LOG_PROGRESS("fit_decorrelation_grid: TPS interpolation requested (Phase 4)");
    /* TODO: Phase 4 - fit TPS surface through grid_phi values */
    /* For now, fall back to bilinear interpolation */
    grid->use_tps_interp = 0;
  }

  LOG_PROGRESS("fit_decorrelation_grid: complete");
  return 0;
}

/**
 * @brief Evaluate decorrelation kernel φ(x, y) at image position.
 *
 * PHASE 3: To be implemented.
 * Uses bilinear or TPS interpolation of grid values.
 */
double eval_decorrelation_at_point(double x, double y,
                                    const decorrelation_grid_t *grid) {
  (void)x;
  (void)y;
  (void)grid;
  LOG_ERROR("eval_decorrelation_at_point: STUB (Phase 3)");
  return -1.0;
}

/**
 * @brief Apply decorrelation to difference image.
 *
 * PHASE 4: To be implemented.
 * Convolves difference image D with spatially-varying φ.
 */
int spatial_convolve_with_decorr(const double *diffImage, int diffImage_ny,
                                  int diffImage_nx, double *diffImageDec,
                                  const decorrelation_grid_t *grid) {
  (void)diffImage;
  (void)diffImage_ny;
  (void)diffImage_nx;
  (void)diffImageDec;
  (void)grid;
  LOG_ERROR("spatial_convolve_with_decorr: STUB (Phase 4)");
  return -1;
}
