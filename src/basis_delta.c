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
 * @brief Delta basis eval_kernel: direct kernel pixel values.
 *
 * @details Reconstructs kernel at (xi, yi) from delta basis coefficients.
 *          Each kernel pixel is directly the fitted coefficient for that basis.
 *          Optional spatial variation (polynomial or TPS) can modulate values.
 *
 * @param[in] xi, yi Image coordinates
 * @param[in] kernelSol Fitted delta basis coefficients (one per kernel pixel)
 * @return Kernel pixel sum
 */
static double delta_eval_kernel(int xi, int yi, double* kernelSol) {
  int basisIdx, kx, ky;
  double coeff, sum;

  if (!kernelSol || fwKernel <= 0 || hwKernel < 0) {
    LOG_ERROR("Invalid parameters for delta kernel evaluation");
    return 0.0;
  }

  sum = 0.0;

  /* Assemble kernel from delta basis coefficients (one per pixel) */
  for (basisIdx = 0; basisIdx < nCompKer; basisIdx++) {
    /* Map basis index to kernel pixel (kx, ky) */
    kx = (basisIdx % fwKernel) - hwKernel;
    ky = (basisIdx / fwKernel) - hwKernel;

    /* Direct coefficient (no spatial variation in simple delta mode) */
    coeff = kernelSol[basisIdx];

    kernel[kx + hwKernel + fwKernel * (ky + hwKernel)] = coeff;
    sum += coeff;
  }

  return sum;
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
    .eval_kernel = delta_eval_kernel,
    .init = delta_init,
    .cleanup = delta_cleanup,
};
