/**
 * @file tps_utils.h
 * @brief Shared thin plate spline (TPS) utility functions for RBF interpolation.
 *
 * Provides core TPS functions used by both Gaussian and delta kernel basis
 * implementations for smooth spatial variation via RBF interpolation.
 *
 * Reference: Bookstein, F.L. "Principal Warps: Thin-Plate Splines and the
 * Decomposition of Deformations". IEEE Trans. PAMI 11(6):567-585, 1989.
 */

#ifndef HOTPANTS_TPS_UTILS_H
#define HOTPANTS_TPS_UTILS_H

#include <math.h>

/**
 * @brief Thin plate spline (TPS) RBF kernel: φ(r) = r² log(r).
 *
 * The TPS kernel is rotation-invariant and minimizes bending energy,
 * making it ideal for smooth spatial interpolation in image differencing.
 *
 * @param r Euclidean distance between two points.
 * @return RBF kernel value φ(r).
 *
 * @note Special case: φ(0) = 0 by convention (r < 1e-10 treated as 0).
 */
static inline double tps_kernel(double r) {
  if (r < 1e-10) return 0.0;
  return r * r * log(r);
}

/**
 * @brief Evaluate a thin plate spline surface at a given point.
 *
 * @details Computes the RBF interpolant:
 *   f(x,y) = c₀ + c₁·x + c₂·y + Σᵢ wᵢ·φ(||pᵢ-(x,y)||)
 *
 * where:
 *   - c₀, c₁, c₂ are polynomial trend coefficients
 *   - pᵢ are control point positions
 *   - wᵢ are RBF weights (fitted coefficients)
 *   - φ is the TPS kernel
 *
 * @param[in] eval_x, eval_y Evaluation point coordinates
 * @param[in] positions Array of n_points (x,y) control point positions (interleaved: x₀,y₀,x₁,y₁,...)
 * @param[in] n_points Number of control points
 * @param[in] weights Array of n_points RBF weights
 * @param[in] poly_coeffs Array of 3 polynomial coefficients [c₀, c₁, c₂]
 * @return Interpolated value at (eval_x, eval_y)
 */
static inline double tps_evaluate(double eval_x, double eval_y, double* positions,
                                   int n_points, double* weights, double* poly_coeffs) {
  int i;
  double value, dx, dy, r;

  /* Initialize with polynomial trend: c₀ + c₁·x + c₂·y */
  value = poly_coeffs[0] + poly_coeffs[1] * eval_x + poly_coeffs[2] * eval_y;

  /* Add RBF contributions: Σᵢ wᵢ·φ(||pᵢ-(x,y)||) */
  for (i = 0; i < n_points; i++) {
    dx = eval_x - positions[2 * i];
    dy = eval_y - positions[2 * i + 1];
    r = sqrt(dx * dx + dy * dy);
    value += weights[i] * tps_kernel(r);
  }

  return value;
}

#endif  /* HOTPANTS_TPS_UTILS_H */
