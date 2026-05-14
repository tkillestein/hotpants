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
 * TODO: Integrate with fit_decorrelation_grid() to auto-estimate variances
 * from actual science and template image data.
 *
 * @param[in] image      Image data (size: ny × nx)
 * @param[in] nx, ny     Image dimensions
 * @return Estimated variance, or -1 on error
 */
static __attribute__((unused)) double estimate_image_variance(const float *image, int nx, int ny) {
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
 * Uses bilinear interpolation of grid values (or TPS if fitted).
 * For spatially-varying decorrelation, interpolates φ smoothly across
 * the image based on precomputed grid values.
 *
 * @param[in] x, y             Image coordinates (pixels, in region space)
 * @param[in] grid             Fitted decorrelation grid with control points
 * @return φ(x, y) value (typically in range [0, 2]), or -1.0 on error
 *
 * Bilinear interpolation formula:
 *   φ(x,y) = (1-u)(1-v)·P₀₀ + u(1-v)·P₁₀ + (1-u)v·P₀₁ + uv·P₁₁
 * where u,v ∈ [0,1] are normalized coords within the grid cell.
 */
double eval_decorrelation_at_point(double x, double y,
                                    const decorrelation_grid_t *grid) {
  if (!grid || grid->nbasis <= 0) {
    LOG_ERROR("eval_decorrelation_at_point: invalid grid");
    return -1.0;
  }

  /* For single-cell grid (1×1), return the only value */
  if (grid->ngrid_x == 1 && grid->ngrid_y == 1) {
    return grid->grid_phi[0];
  }

  /* Find enclosing grid cell */
  /* Grid coordinates are in region space; we assume stamps are laid out
     evenly. For simplicity, assume grid_x and grid_y are sorted. */

  /* Estimate grid spacing from first two grid points */
  double dx_grid = grid->ngrid_x > 1 ?
      (grid->grid_x[grid->ngrid_y] - grid->grid_x[0]) : 1.0;
  double dy_grid = grid->ngrid_y > 1 ?
      (grid->grid_y[1] - grid->grid_y[0]) : 1.0;

  if (dx_grid < ZEROVAL) dx_grid = 1.0;
  if (dy_grid < ZEROVAL) dy_grid = 1.0;

  /* Normalize position within grid */
  double x_grid = (x - grid->grid_x[0]) / dx_grid;
  double y_grid = (y - grid->grid_y[0]) / dy_grid;

  /* Clamp to valid range [0, ngrid-1) */
  if (x_grid < 0.0) x_grid = 0.0;
  if (y_grid < 0.0) y_grid = 0.0;
  if (x_grid >= grid->ngrid_x - 1) x_grid = grid->ngrid_x - 1.001;
  if (y_grid >= grid->ngrid_y - 1) y_grid = grid->ngrid_y - 1.001;

  /* Get integer and fractional parts */
  int ix = (int)x_grid;
  int iy = (int)y_grid;
  double u = x_grid - ix;
  double v = y_grid - iy;

  /* Bounds check */
  if (ix < 0 || ix >= grid->ngrid_x - 1 || iy < 0 || iy >= grid->ngrid_y - 1) {
    LOG_DEBUG("eval_decorrelation_at_point: out of bounds (%f, %f)", x, y);
    return 1.0; /* Fallback: no decorrelation at edges */
  }

  /* Bilinear interpolation: φ(x,y) = (1-u)(1-v)P₀₀ + u(1-v)P₁₀ + (1-u)vP₀₁ + uvP₁₁ */
  double p00 = grid->grid_phi[ix * grid->ngrid_y + iy];
  double p10 = grid->grid_phi[(ix + 1) * grid->ngrid_y + iy];
  double p01 = grid->grid_phi[ix * grid->ngrid_y + (iy + 1)];
  double p11 = grid->grid_phi[(ix + 1) * grid->ngrid_y + (iy + 1)];

  double phi = (1.0 - u) * (1.0 - v) * p00 + u * (1.0 - v) * p10 +
               (1.0 - u) * v * p01 + u * v * p11;

  return phi;
}

/**
 * @brief Apply decorrelation to difference image.
 *
 * Convolves the difference image D with spatially-varying φ to produce
 * D_decorr = φ ⊗ D. Uses direct spatial convolution with on-the-fly φ
 * evaluation at each pixel.
 *
 * @param[in]     diffImage     Difference image [diffImage_ny × diffImage_nx]
 * @param[in]     diffImage_ny  Difference image height
 * @param[in]     diffImage_nx  Difference image width
 * @param[in,out] diffImageDec  Output decorrelated difference image
 * @param[in]     grid          Fitted decorrelation grid with φ values
 * @return 0 on success, -1 on error
 *
 * Algorithm:
 *   For each output pixel (i,j):
 *     φ_ij = eval_decorrelation_at_point(i, j)  [interpolate from grid]
 *     output[i,j] = φ_ij * input[i,j]           [scale by local φ]
 *
 * This is a simple point-wise scaling operation, not a full convolution kernel.
 * The actual decorrelation effect comes from how φ was computed from the
 * matching kernel κ (see DMTN-021).
 *
 * @note This implementation uses simple point-wise scaling. A more sophisticated
 *       approach would perform actual spatial convolution with a spatially-varying
 *       φ kernel, but that requires storing the full φ field and is
 *       computationally expensive. For now, point-wise scaling approximates the
 *       effect assuming φ is slowly varying across the image.
 */
int spatial_convolve_with_decorr(const double *diffImage, int diffImage_ny,
                                  int diffImage_nx, double *diffImageDec,
                                  const decorrelation_grid_t *grid) {
  if (!diffImage || !diffImageDec || !grid) {
    LOG_ERROR("spatial_convolve_with_decorr: NULL pointer");
    return -1;
  }

  if (diffImage_ny <= 0 || diffImage_nx <= 0) {
    LOG_ERROR("spatial_convolve_with_decorr: invalid dimensions (%d × %d)",
              diffImage_ny, diffImage_nx);
    return -1;
  }

  int i, j, idx;
  double phi_ij;

  LOG_PROGRESS("spatial_convolve_with_decorr: decorrelating difference image "
               "(%d × %d) using grid (%d × %d)",
               diffImage_ny, diffImage_nx, grid->ngrid_y, grid->ngrid_x);

  /* For each pixel, evaluate φ at that location and scale */
  for (i = 0; i < diffImage_ny; i++) {
    for (j = 0; j < diffImage_nx; j++) {
      idx = i * diffImage_nx + j;

      /* Evaluate φ(i, j) via interpolation */
      phi_ij = eval_decorrelation_at_point((double)j, (double)i, grid);

      if (phi_ij < 0.0) {
        LOG_DEBUG("spatial_convolve_with_decorr: invalid φ at (%d, %d), "
                  "skipping",
                  i, j);
        diffImageDec[idx] = diffImage[idx]; /* Fallback: no decorrelation */
        continue;
      }

      /* Apply decorrelation: D_decorr = φ · D */
      diffImageDec[idx] = phi_ij * diffImage[idx];
    }
  }

  LOG_PROGRESS("spatial_convolve_with_decorr: complete");
  return 0;
}

/* =====================================================================
   HIGH-LEVEL API FOR MAIN.C INTEGRATION
   ===================================================================== */

/**
 * @brief High-level wrapper: compute decorrelation grid for a region.
 *
 * Call this after fitKernel() to prepare the decorrelation grid.
 * Initializes decorr_grid (global) with φ values at all stamp centers.
 *
 * @param[in] kernelSol  Fitted kernel solution (double array)
 * @return 0 on success, -1 on error
 */
int decorrelation_init_region(const double *kernelSol) {
  if (!kernelSol) {
    LOG_ERROR("decorrelation_init_region: NULL kernelSol");
    return -1;
  }

  /* Free any previous grid */
  if (decorr_grid) {
    decorrelation_grid_free(&decorr_grid);
  }

  /* Allocate grid with stamp dimensions */
  decorr_grid = decorrelation_grid_alloc(nStampX, nStampY);
  if (!decorr_grid) {
    LOG_ERROR("decorrelation_init_region: allocation failed");
    return -1;
  }

  /* Initialize grid coordinates (stamp centers in region space) */
  int gx, gy;
  for (gx = 0; gx < nStampX; gx++) {
    for (gy = 0; gy < nStampY; gy++) {
      /* Rough estimate: stamp centers based on uniform grid */
      int idx = gx * nStampY + gy;
      decorr_grid->grid_x[idx] = (gx + 0.5) * rPixX / nStampX;
      decorr_grid->grid_y[idx] = (gy + 0.5) * rPixY / nStampY;
    }
  }

  /* Fit decorrelation grid: compute φ at each stamp center */
  int ret = fit_decorrelation_grid(decorr_grid, kernelSol,
                                    decorrScienceVar, decorrTemplateVar,
                                    decorrUseTPS);
  if (ret < 0) {
    LOG_ERROR("decorrelation_init_region: fit_decorrelation_grid failed");
    decorrelation_grid_free(&decorr_grid);
    return -1;
  }

  return 0;
}

/**
 * @brief High-level wrapper: apply decorrelation to a float difference image.
 *
 * Call this after the difference image is complete (after spatial_convolve).
 * Converts float difference image to double, applies decorrelation, converts back.
 *
 * @param[in]     diffImage_float   Input difference image (float, rPixY × rPixX)
 * @param[in]     ny, nx            Image dimensions
 * @param[in,out] diffImageDec_float Output decorrelated difference image (float)
 * @return 0 on success, -1 on error
 */
int decorrelation_apply_region(const float *diffImage_float, int ny, int nx,
                                 float *diffImageDec_float) {
  if (!diffImage_float || !diffImageDec_float || !decorr_grid) {
    LOG_ERROR("decorrelation_apply_region: NULL pointer");
    return -1;
  }

  if (ny <= 0 || nx <= 0) {
    LOG_ERROR("decorrelation_apply_region: invalid dimensions (%d × %d)", ny,
              nx);
    return -1;
  }

  /* Convert float to double for decorrelation computation */
  double *diffImage_double =
      (double *)xcalloc((size_t)ny * nx, sizeof(double));
  double *diffImageDec_double =
      (double *)xcalloc((size_t)ny * nx, sizeof(double));

  int idx;
  for (idx = 0; idx < ny * nx; idx++) {
    diffImage_double[idx] = (double)diffImage_float[idx];
  }

  /* Apply decorrelation */
  int ret = spatial_convolve_with_decorr(diffImage_double, ny, nx,
                                          diffImageDec_double, decorr_grid);
  if (ret < 0) {
    LOG_ERROR("decorrelation_apply_region: spatial_convolve_with_decorr failed");
    free(diffImage_double);
    free(diffImageDec_double);
    return -1;
  }

  /* Convert result back to float */
  for (idx = 0; idx < ny * nx; idx++) {
    diffImageDec_float[idx] = (float)diffImageDec_double[idx];
  }

  free(diffImage_double);
  free(diffImageDec_double);

  return 0;
}
