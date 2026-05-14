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
 * PHASE 2: To be implemented.
 * Reconstructs K(gx, gy), FFTs it, applies decorrelation formula.
 */
int compute_phi_at_gridpoint(int gx, int gy, const double *kernelSol,
                              decorrelation_grid_t *grid, double sig_s2,
                              double sig_t2) {
  (void)gx;
  (void)gy;
  (void)kernelSol;
  (void)grid;
  (void)sig_s2;
  (void)sig_t2;
  LOG_ERROR("compute_phi_at_gridpoint: STUB (Phase 2)");
  return -1;
}

/**
 * @brief Fit decorrelation control grid.
 *
 * PHASE 2: To be implemented.
 * Iterates over stamp grid, reconstructs kernel at each point, applies formula.
 */
int fit_decorrelation_grid(decorrelation_grid_t *grid,
                            const double *kernelSol, double sig_s2,
                            double sig_t2, int use_tps) {
  (void)grid;
  (void)kernelSol;
  (void)sig_s2;
  (void)sig_t2;
  (void)use_tps;
  LOG_ERROR("fit_decorrelation_grid: STUB (Phase 2)");
  return -1;
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
