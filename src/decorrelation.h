/**
 * @file decorrelation.h
 * @brief Afterburner decorrelation for Alard & Lupton image subtraction.
 *
 * Implements the LSST DMTN-021 decorrelation algorithm to remove correlated
 * noise introduced by PSF-matching convolution.
 *
 * For non-spatially-varying kernels:
 *   φ̂(k) = √[(σ_s² + σ_t²) / (σ_s² + σ_t² · κ̂²(k))]
 *
 * For spatially-varying kernels (polynomial or TPS), φ is evaluated on a control
 * grid (stamp centers) and interpolated spatially. The difference image is then
 * convolved with φ to produce D_decorr = φ ⊗ D.
 */

#ifndef HOTPANTS_DECORRELATION_H
#define HOTPANTS_DECORRELATION_H

#include <fftw3.h>

/* =====================================================================
   DECORRELATION GRID DATA STRUCTURE
   ===================================================================== */

/**
 * @brief Control grid for spatially-varying decorrelation kernel.
 *
 * Stores precomputed φ values at stamp centers and supports spatial
 * interpolation to image pixels. For single-region mode, the grid matches
 * the stamp grid dimensions (nStampX × nStampY).
 */
typedef struct decorrelation_grid {
    int ngrid_x, ngrid_y;               /**< Grid dimensions (matches stamp grid) */
    int nbasis;                         /**< Total grid points (ngrid_x * ngrid_y) */

    double *grid_phi;                   /**< φ values at grid points [nbasis] */
    double *grid_x, *grid_y;            /**< Grid center coordinates [nbasis] */

    /* Optional: TPS interpolation of φ (for smoother spatial variation) */
    double *phi_tps_weights;            /**< RBF weights for TPS [nbasis] */
    double *phi_tps_poly;               /**< Polynomial trend [3] */
    int use_tps_interp;                 /**< 1 = use TPS, 0 = bilinear interp */

    /* Variance estimates */
    double sigma_science2;              /**< σ_s² (science image variance) */
    double sigma_template2;             /**< σ_t² (template image variance) */

    /* FFT workspace (cached for multiple convolutions) */
    fftw_plan fft_plan;                 /**< FFTW plan for φ convolution */
    double *phi_fft;                    /**< FFT of spatially-averaged φ [fwKernel²] */
} decorrelation_grid_t;

/* =====================================================================
   FUNCTION SIGNATURES
   ===================================================================== */

/**
 * @brief Compute decorrelation kernel φ̂(k) in Fourier space.
 *
 * Applies the formula:
 *   φ̂(k) = √[(σ_s² + σ_t²) / (σ_s² + σ_t² · κ̂²(k))]
 *
 * @param[in]  kappa_ft        PSF-matching kernel in Fourier space [fwKernel²]
 * @param[in]  sig_s2          Science image variance σ_s²
 * @param[in]  sig_t2          Template image variance σ_t²
 * @param[out] phi_ft          Output decorrelation kernel [fwKernel²]
 * @param[in]  fwKernel        Kernel full-width (assumes square kernel)
 * @return 0 on success, -1 on numerical error
 *
 * @details
 * Handles numerical edge cases:
 * - Avoids division by zero if denominator is very small
 * - Ensures φ ≥ 0 (takes absolute value if needed)
 * - Clamps extreme values to prevent overflow
 */
int compute_decorrelation_kernel_fft(const double *kappa_ft, double sig_s2,
                                      double sig_t2, double *phi_ft, int fwKernel);

/**
 * @brief Reconstruct local kernel and compute its decorrelation kernel.
 *
 * At grid point (gx, gy), reconstructs the spatially-varying kernel K(gx, gy),
 * FFTs it, applies the decorrelation formula, and stores the result.
 *
 * @param[in]  gx, gy          Grid point index (0 to ngrid_x-1, 0 to ngrid_y-1)
 * @param[in]  kernelSol       Fitted kernel coefficients from fitKernel()
 * @param[in]  grid            Decorrelation grid (output stored in grid->grid_phi[gx*ngrid_y+gy])
 * @param[in]  sig_s2          Science image variance
 * @param[in]  sig_t2          Template image variance
 * @return 0 on success, -1 on error
 */
int compute_phi_at_gridpoint(int gx, int gy, const double *kernelSol,
                              decorrelation_grid_t *grid,
                              double sig_s2, double sig_t2);

/**
 * @brief Fit decorrelation control grid.
 *
 * Pre-computes φ values at all stamp centers by reconstructing the local kernel
 * and applying the decorrelation formula. Optionally fits a TPS surface through
 * the grid points for smooth spatial interpolation.
 *
 * @param[in,out] grid          Allocated decorrelation grid
 * @param[in]     kernelSol     Fitted kernel coefficients from fitKernel()
 * @param[in]     sig_s2        Science image variance (or 0 to auto-estimate)
 * @param[in]     sig_t2        Template image variance (or 0 to auto-estimate)
 * @param[in]     use_tps       1 = fit TPS, 0 = use bilinear interpolation
 * @return 0 on success, -1 on error
 *
 * @details
 * If sig_s2 or sig_t2 are 0, estimates them from image statistics.
 * Fits optional TPS surface through grid_phi values for smooth evaluation.
 */
int fit_decorrelation_grid(decorrelation_grid_t *grid, const double *kernelSol,
                            double sig_s2, double sig_t2, int use_tps);

/**
 * @brief Evaluate decorrelation kernel φ(x, y) at image position.
 *
 * Uses bilinear or TPS interpolation of grid values to evaluate φ at arbitrary
 * pixel position (x, y) in region coordinates.
 *
 * @param[in] x, y             Image coordinates (pixels)
 * @param[in] grid             Fitted decorrelation grid
 * @return φ(x, y) value, or -1.0 if interpolation fails
 */
double eval_decorrelation_at_point(double x, double y,
                                    const decorrelation_grid_t *grid);

/**
 * @brief Apply decorrelation to difference image.
 *
 * Convolves the difference image D with spatially-varying φ to produce
 * D_decorr = φ ⊗ D. Uses FFT-accelerated convolution.
 *
 * @param[in]     diffImage     Difference image [rPixY × rPixX]
 * @param[in]     diffImage_ny  Difference image height
 * @param[in]     diffImage_nx  Difference image width
 * @param[in,out] diffImageDec  Output decorrelated difference image [rPixY × rPixX]
 * @param[in]     grid          Fitted decorrelation grid
 * @return 0 on success, -1 on error
 *
 * @details
 * For spatially-varying φ, uses direct spatial convolution with on-the-fly
 * φ evaluation. For non-varying kernels, uses cached FFT.
 */
int spatial_convolve_with_decorr(const double *diffImage, int diffImage_ny,
                                  int diffImage_nx, double *diffImageDec,
                                  const decorrelation_grid_t *grid);

/**
 * @brief Allocate decorrelation grid.
 *
 * @param[in] ngrid_x, ngrid_y Grid dimensions
 * @return Pointer to allocated grid, or NULL on error
 */
decorrelation_grid_t* decorrelation_grid_alloc(int ngrid_x, int ngrid_y);

/**
 * @brief Free decorrelation grid and all associated memory.
 *
 * @param[in,out] grid  Grid to free (pointer set to NULL on return)
 * @return 0 on success
 */
int decorrelation_grid_free(decorrelation_grid_t **grid);

#endif  /* HOTPANTS_DECORRELATION_H */
