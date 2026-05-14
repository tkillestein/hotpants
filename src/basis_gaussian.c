/**
 * @file basis_gaussian.c
 * @brief Gaussian kernel basis implementation (multi-Gaussian with polynomial variation).
 *
 * References:
 * - Alard & Lupton (1998), "A Method for Optimal Image Subtraction", ApJ 503:325
 * - Section 2.2: Kernel parameterization as sum of fixed Gaussian basis functions
 *   weighted by spatial polynomials
 */

#include <string.h>
#include <math.h>
#include "defaults.h"
#include "globals.h"
#include "basis.h"
#include "functions.h"

/* =====================================================================
   GAUSSIAN BASIS IMPLEMENTATION
   ===================================================================== */

/**
 * @brief Gaussian basis convolve_stamp: correlate with precomputed basis.
 *
 * @details Uses precomputed kernel_vec[basisIdx] which is a fixed Gaussian basis
 *          function. Performs separable 2D convolution via xy_conv_stamp().
 *
 * @param[in] stamp Stamp array (mStampX × mStampY)
 * @param[in] mStampX, mStampY Stamp dimensions
 * @param[in] basisIdx Basis function index (0 to nbasis-1)
 * @param[out] pixFitted Output correlation array (mStampX × mStampY)
 * @return 0 on success, -1 on error
 */
static int gaussian_convolve_stamp(double* stamp, int mStampX, int mStampY,
                                   int basisIdx, double* pixFitted) {
  (void)stamp;      /* Unused: Gaussian uses precomputed kernel_vec */
  (void)mStampX;
  (void)mStampY;
  (void)pixFitted;

  if (basisIdx < 0 || basisIdx >= nCompKer) {
    LOG_ERROR("Invalid basisIdx %d (nCompKer=%d)", basisIdx, nCompKer);
    return -1;
  }

  /* For Gaussian basis, convolution is already computed in stamp->vectors[]
     during fillStamp() via xy_conv_stamp(). This function is a no-op for
     Gaussian basis since vectors are pre-filled. */
  LOG_DEBUG("Gaussian basis: using precomputed stamp->vectors[%d]", basisIdx);

  return 0;
}

/**
 * @brief Gaussian basis eval_kernel: polynomial-weighted Gaussian sum.
 *
 * @details Evaluates spatial polynomial weights at (xi, yi) and assembles
 *          the kernel as: K(x,y) = Σᵢ cᵢ(x,y) · φᵢ(x,y)
 *          where cᵢ are polynomial-weighted coefficients and φᵢ are Gaussians.
 *
 * @param[in] xi, yi Image coordinates
 * @param[in] kernelSol Fitted Gaussian coefficients
 * @return Kernel pixel sum
 */
static double gaussian_eval_kernel(int xi, int yi, double* kernelSol) {
  int gaussianCompIdx, solutionIdx, polyDegX, polyDegY, kernelPixelIdx;
  int pixelCompIdx;
  double polyBasisX, polyBasisY, kernelSum;
  double normalizedX, normalizedY;

  solutionIdx = 2;
  /* RANGE FROM -1 to 1 */
  normalizedX = (xi - 0.5 * rPixX) / (0.5 * rPixX);
  normalizedY = (yi - 0.5 * rPixY) / (0.5 * rPixY);

  /* Compute spatially-varying kernel coefficients from polynomials */
  for (gaussianCompIdx = 1; gaussianCompIdx < nCompKer; gaussianCompIdx++) {
    kernel_coeffs[gaussianCompIdx] = 0.0;
    polyBasisX = 1.0;
    for (polyDegX = 0; polyDegX <= kerOrder; polyDegX++) {
      polyBasisY = 1.0;
      for (polyDegY = 0; polyDegY <= kerOrder - polyDegX; polyDegY++) {
        kernel_coeffs[gaussianCompIdx] +=
            kernelSol[solutionIdx++] * polyBasisX * polyBasisY;
        polyBasisY *= normalizedY;
      }
      polyBasisX *= normalizedX;
    }
  }
  kernel_coeffs[0] = kernelSol[1];

  /* Assemble kernel from weighted Gaussian basis functions */
  for (kernelPixelIdx = 0; kernelPixelIdx < fwKernel * fwKernel;
       kernelPixelIdx++)
    kernel[kernelPixelIdx] = 0.0;

  kernelSum = 0.0;
  for (kernelPixelIdx = 0; kernelPixelIdx < fwKernel * fwKernel;
       kernelPixelIdx++) {
    for (pixelCompIdx = 0; pixelCompIdx < nCompKer; pixelCompIdx++) {
      kernel[kernelPixelIdx] += kernel_coeffs[pixelCompIdx] *
                                kernel_vec[pixelCompIdx][kernelPixelIdx];
    }
    kernelSum += kernel[kernelPixelIdx];
  }
  return kernelSum;
}

/**
 * @brief Gaussian basis init: set up Gaussian basis functions.
 *
 * @details Calls getKernelVec() to create multi-Gaussian basis.
 *          Sets nCompKer to total number of basis functions.
 *
 * @return nCompKer on success, -1 on error
 */
static int gaussian_init(void) {
  int nbasis;

  LOG_PROGRESS("Initializing Gaussian basis...");

  nbasis = getKernelVec();
  if (nbasis < 0) {
    LOG_ERROR("Failed to create Gaussian basis functions");
    return -1;
  }

  nCompKer = nbasis;
  LOG_PROGRESS("Gaussian basis: %d basis functions", nCompKer);

  return nCompKer;
}

/**
 * @brief Gaussian basis cleanup: free basis-specific resources.
 *
 * @details Frees kernel_vec arrays. Called on basis switch or exit.
 *
 * @return 0 on success
 */
static int gaussian_cleanup(void) {
  LOG_DEBUG("Cleaning up Gaussian basis");

  /* kernel_vec is freed in main cleanup path; no additional cleanup needed */
  return 0;
}

/* =====================================================================
   GAUSSIAN BASIS REGISTRATION
   ===================================================================== */

kernel_basis_t gaussian_basis = {
    .name = "gaussian",
    .nbasis = 0,  /* Set by init() */
    .convolve_stamp = gaussian_convolve_stamp,
    .eval_kernel = gaussian_eval_kernel,
    .init = gaussian_init,
    .cleanup = gaussian_cleanup,
};
