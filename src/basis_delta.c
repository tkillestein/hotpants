/**
 * @file basis_delta.c
 * @brief Delta function kernel basis implementation (pixel-level delta functions).
 *
 * References:
 * - Bramich (2008), "The Optimal Difference Image Combination of Dithered
 *   Images", MNRAS 389:1365 (arXiv:0802.1273)
 * - Section 2.1: Delta function basis with Laplacian regularization
 *
 * Each kernel pixel becomes a separate basis function. For fwKernel×fwKernel
 * kernel, there are fwKernel² basis functions (one per pixel).
 */

#include <string.h>
#include <math.h>
#include "defaults.h"
#include "globals.h"
#include "basis.h"
#include "functions.h"
#include "tps_utils.h"

/* =====================================================================
   DELTA BASIS IMPLEMENTATION
   ===================================================================== */

/**
 * @brief Delta basis convolve_stamp: extract stamp value at delta position.
 *
 * @details For delta basis function at kernel pixel (kx, ky), extracts the
 *          stamp value at each position (stamp_x, stamp_y), zero-padded
 *          outside bounds.
 *
 * @param[in] stamp Stamp pixel array (mStampX × mStampY)
 * @param[in] mStampX, mStampY Stamp dimensions
 * @param[in] basisIdx Kernel pixel index (0 to nbasis-1)
 * @param[out] pixFitted Output array (mStampX × mStampY)
 * @return 0 on success, -1 on error
 */
static int delta_convolve_stamp(double* stamp, int mStampX, int mStampY,
                                int basisIdx, double* pixFitted) {
  int stampPixelX, stampPixelY, pixelIdx;
  int kernel_x, kernel_y;

  if (basisIdx < 0 || basisIdx >= fwKernel * fwKernel) {
    LOG_ERROR("Invalid basisIdx %d (kernel has %d×%d = %d pixels)", basisIdx,
              fwKernel, fwKernel, fwKernel * fwKernel);
    return -1;
  }

  /* Map basis index to kernel offset (kx, ky) */
  kernel_x = (basisIdx % fwKernel) - hwKernel;
  kernel_y = (basisIdx / fwKernel) - hwKernel;

  /* Extract stamp values shifted by kernel offset (zero-padded outside bounds) */
  pixelIdx = 0;
  for (stampPixelY = 0; stampPixelY < mStampY; stampPixelY++) {
    for (stampPixelX = 0; stampPixelX < mStampX; stampPixelX++) {
      int stamp_y = stampPixelY + kernel_y;
      int stamp_x = stampPixelX + kernel_x;

      /* Zero-pad outside stamp bounds */
      if (stamp_x >= 0 && stamp_x < mStampX && stamp_y >= 0 &&
          stamp_y < mStampY) {
        int stamp_idx = stamp_x + mStampX * stamp_y;
        pixFitted[pixelIdx] = stamp[stamp_idx];
      } else {
        pixFitted[pixelIdx] = 0.0;
      }
      pixelIdx++;
    }
  }

  return 0;
}

/**
 * @brief Delta basis eval_kernel: kernel pixel values with optional spatial variation.
 *
 * @details Reconstructs kernel at (xi, yi) from delta basis coefficients.
 *          For delta basis, kernelSol stores polynomial coefficients for each
 *          kernel pixel (one coefficient per pixel), with spatial variation applied
 *          via polynomial weighting:
 *
 *            K_pixel(xi, yi) = Σ_{ix,iy} kernelSol[pixel_idx][ix,iy] · x^ix · y^iy
 *
 *          For non-spatial mode (kerOrder=0), each pixel has a single coefficient.
 *
 * @param[in] xi, yi Image coordinates
 * @param[in] kernelSol Fitted delta basis coefficients (with polynomial spatial variation)
 * @return Kernel pixel sum
 */
static double delta_eval_kernel(int xi, int yi, double* kernelSol) {
  int basisIdx, kx, ky, kernelPixelIdx, solutionIdx;
  double coeff, sum, normalizedX, normalizedY, polyBasisX, polyBasisY;
  int polyDegX, polyDegY;

  if (!kernelSol || fwKernel <= 0 || hwKernel < 0) {
    LOG_ERROR("Invalid parameters for delta kernel evaluation");
    return 0.0;
  }

  sum = 0.0;

  /* Normalize coordinates to [-1, 1] for polynomial evaluation */
  normalizedX = (xi - 0.5 * rPixX) / (0.5 * rPixX);
  normalizedY = (yi - 0.5 * rPixY) / (0.5 * rPixY);

  /* For delta basis with spatial variation:
     kernelSol layout is [const, coeff_1_0, coeff_1_1, ..., coeff_N_polys...]
     where each pixel has (kerOrder+1)*(kerOrder+2)/2 polynomial terms.

     Note: kernelSol[0] contains initial background/flux scaling (not used for kernel pixels).
  */

  solutionIdx = 1;  /* Skip initial term */

  /* Assemble kernel from delta basis coefficients with polynomial spatial variation */
  for (basisIdx = 0; basisIdx < nCompKer; basisIdx++) {
    /* Map basis index to kernel pixel (kx, ky) */
    kx = (basisIdx % fwKernel) - hwKernel;
    ky = (basisIdx / fwKernel) - hwKernel;
    kernelPixelIdx = kx + hwKernel + fwKernel * (ky + hwKernel);

    /* Evaluate polynomial spatial variation for this kernel pixel */
    coeff = 0.0;
    polyBasisX = 1.0;
    for (polyDegX = 0; polyDegX <= kerOrder; polyDegX++) {
      polyBasisY = 1.0;
      for (polyDegY = 0; polyDegY <= kerOrder - polyDegX; polyDegY++) {
        if (solutionIdx >= 0) {  /* Bounds check */
          coeff += kernelSol[solutionIdx] * polyBasisX * polyBasisY;
        }
        solutionIdx++;
        polyBasisY *= normalizedY;
      }
      polyBasisX *= normalizedX;
    }

    kernel[kernelPixelIdx] = coeff;
    sum += coeff;
  }

  return sum;
}

/* tps_kernel() and tps_evaluate() are now in tps_utils.h (included above) */

/**
 * @brief Get offset in kernelSol for polynomial (non-TPS) coefficients.
 *
 * Layout (same for all basis types):
 *   [0]: reserved
 *   [1..nComp]: kernel polynomial coefficients
 *   [nComp+1..nComp+nbg_vec]: background polynomial coefficients
 *
 * where nComp = (nCompKer-1) * (kerOrder+1)*(kerOrder+2)/2
 *       nbg_vec = (bgOrder+1)*(bgOrder+2)/2
 */
static int delta_kernelSol_size_polynomial(void) {
  /* All nCompKer pixels have polynomial variation; no DC slot.
   * poly_size = mat_size + 1 = nCompKer*ncomp2 + nbg_vec + 1. */
  int ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  int ncomp = nCompKer * ncomp2;
  int nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
  return ncomp + nbg_vec + 1;
}

/**
 * @brief Get offset in kernelSol for TPS RBF weights of a kernel component.
 *
 * For delta basis, component indices are kernel pixels (0..nCompKer-1).
 *
 * @param comp_idx Kernel component index (0..nCompKer-1)
 * @param n_stamps Number of stamps (control points)
 * @return Offset in kernelSol; access weights via kernelSol[offset..offset+n_stamps)
 */
static int delta_kernelSol_offset_tps_weights(int comp_idx, int n_stamps) {
  int poly_size = delta_kernelSol_size_polynomial();
  return poly_size + comp_idx * n_stamps;
}

/**
 * @brief Get offset in kernelSol for TPS polynomial trend of a kernel component.
 *
 * @param comp_idx Kernel component index (0..nCompKer-1)
 * @param n_stamps Number of stamps
 * @return Offset in kernelSol; access 3 poly coeffs via kernelSol[offset..offset+3)
 */
static int delta_kernelSol_offset_tps_poly(int comp_idx, int n_stamps) {
  int poly_size = delta_kernelSol_size_polynomial();
  return poly_size + nCompKer * n_stamps + comp_idx * 3;
}

/**
 * @brief Get offset in kernelSol for stamp center positions.
 *
 * @param n_stamps Number of stamps
 * @return Offset in kernelSol; positions are stored as x0,y0,x1,y1,...
 */
static int delta_kernelSol_offset_tps_positions(int n_stamps) {
  int poly_size = delta_kernelSol_size_polynomial();
  return poly_size + nCompKer * n_stamps + 3 * nCompKer;
}

/**
 * @brief Delta basis eval_kernel with TPS spatial variation.
 *
 * @details For delta basis with TPS, each kernel pixel coefficient is represented
 *          as an RBF interpolant over the stamp centers:
 *
 *            K_pixel(xi, yi) = c₀ + c₁·xi + c₂·yi + Σⱼ wⱼ·φ(||stampⱼ-(xi,yi)||)
 *
 *          where:
 *          - (xi, yi) is the evaluation point (image position)
 *          - stampⱼ are the fitted stamp center positions (control points)
 *          - wⱼ are RBF weights for this pixel
 *          - c₀, c₁, c₂ are polynomial trend coefficients
 *          - φ is the thin plate spline kernel φ(r) = r² log(r)
 *
 *          This provides smooth spatial variation of each kernel pixel value
 *          without tiling artifacts (suitable for single-region processing).
 *
 * @param[in] xi, yi Image coordinates where kernel is evaluated
 * @param[in] kernelSol Extended kernel solution vector (includes TPS parameters)
 * @return Sum of all kernel pixels (flux scaling factor)
 */
static double delta_eval_kernel_tps(int xi, int yi, double* kernelSol) {
  int basisIdx, kx, ky, kernelPixelIdx;
  double coeff, sum, *tps_weights, *tps_poly, *positions;
  int pos_offset;

  if (!kernelSol || fwKernel <= 0 || hwKernel < 0) {
    LOG_ERROR("Invalid parameters for delta+TPS kernel evaluation");
    return 0.0;
  }

  sum = 0.0;

  /* Get stamp center positions (control points for RBF interpolation) */
  pos_offset = delta_kernelSol_offset_tps_positions(nS);
  positions = kernelSol + pos_offset;

  LOG_DEBUG("Delta+TPS: evaluating kernel with %d pixels, %d control points", nCompKer, nS);

  /* For each kernel pixel (basis function), evaluate its TPS-interpolated value
     at the image position (xi, yi). This gives smooth spatial variation of each
     pixel without polynomial discontinuities. */
  for (basisIdx = 0; basisIdx < nCompKer; basisIdx++) {
    /* Map basis index to kernel pixel (kx, ky) */
    kx = (basisIdx % fwKernel) - hwKernel;
    ky = (basisIdx / fwKernel) - hwKernel;
    kernelPixelIdx = kx + hwKernel + fwKernel * (ky + hwKernel);

    /* Get TPS weights and polynomial coefficients for this pixel */
    tps_weights = kernelSol + delta_kernelSol_offset_tps_weights(basisIdx, nS);
    tps_poly = kernelSol + delta_kernelSol_offset_tps_poly(basisIdx, nS);

    /* Evaluate RBF surface: value = c₀ + c₁·xi + c₂·yi + Σⱼ wⱼ·φ(||posⱼ-(xi,yi)||) */
    coeff = tps_evaluate((double)xi, (double)yi, positions, nS, tps_weights, tps_poly);

    kernel[kernelPixelIdx] = coeff;
    sum += coeff;
  }

  LOG_DEBUG("Delta+TPS: kernel sum = %.6f", sum);
  return sum;
}

/**
 * @brief Dispatch delta kernel evaluation to polynomial or TPS variant.
 *
 * @param[in] xi, yi Image coordinates
 * @param[in] kernelSol Kernel solution (layout depends on useTPS)
 * @return Kernel pixel sum
 */
double delta_eval_kernel_dispatch(int xi, int yi, double* kernelSol) {
  if (useTPS) {
    return delta_eval_kernel_tps(xi, yi, kernelSol);
  } else {
    return delta_eval_kernel(xi, yi, kernelSol);
  }
}

/**
 * @brief Delta basis init: set up pixel-level delta functions.
 *
 * @details Initializes nCompKer = fwKernel² (one basis function per kernel pixel).
 *          No dynamic allocation needed; basis functions are implicit.
 *
 * @return nCompKer on success, -1 on error
 */
static int delta_init(void) {
  if (hwKernel <= 0 || fwKernel <= 0) {
    LOG_ERROR("hwKernel must be positive; cannot initialize delta basis");
    return -1;
  }

  nCompKer = fwKernel * fwKernel;

  LOG_PROGRESS("Delta function basis: %d×%d = %d basis functions (one per kernel pixel)",
               fwKernel, fwKernel, nCompKer);

  return nCompKer;
}

/**
 * @brief Delta basis cleanup: free basis-specific resources.
 *
 * @details Delta basis requires no dynamic allocation, so cleanup is empty.
 *
 * @return 0 on success
 */
static int delta_cleanup(void) {
  LOG_DEBUG("Cleaning up Delta basis");
  return 0;
}

/* =====================================================================
   DELTA BASIS REGISTRATION
   ===================================================================== */

kernel_basis_t delta_basis = {
    .name = "delta",
    .nbasis = 0,  /* Set by init() */
    .convolve_stamp = delta_convolve_stamp,
    .eval_kernel = delta_eval_kernel_dispatch,
    .init = delta_init,
    .cleanup = delta_cleanup,
};
