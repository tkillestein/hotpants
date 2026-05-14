#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fitsio.h>
#include <lapacke.h>
#include <fftw3.h>

#include "defaults.h"
#include "globals.h"
#include "functions.h"

/*

Several of these subroutines appear originally in code created by
Cristophe Alard for ISIS, but have been modified and/or rewritten
for the current software package.  In particular, the construction
of the least squares matrices have been taken directly from the ISIS
code.

08/20/01 acbecker@physics.bell-labs.com

*/

/* Forward declarations for kernel evaluation dispatchers and TPS functions */
double make_kernel_tps(int xi, int yi, double* kernelSol);
double make_kernel_dispatch(int xi, int yi, double* kernelSol);
static double make_kernel_local_dispatch(int xi, int yi, double* kernelSol,
                                         double* lkernel, double* lkernel_coeffs);
static int tps_fit_background(stamp_struct* stamps, int n_stamps,
                              double* kernelSol);
double get_background_tps(int xi, int yi, double* kernelSol);

/* Forward declarations for delta function basis (Phase 2+) */
/* TODO: Implement delta basis functions in Phase 2
double make_kernel_delta(int xi, int yi, double* kernelSol);
int xy_conv_stamp_delta(double* stamp, int mStampX, int mStampY, int basisIdx,
                        double* pixFitted);
int build_matrix_delta(int substampIdx, double* kernelSol, int* iMatrixSize);
*/

/* =====================================================================
   THIN PLATE SPLINE (TPS) SPATIAL VARIATION — Core RBF Functions
   ===================================================================== */

/**
 * @brief Thin plate spline (TPS) RBF kernel: φ(r) = r² log(r).
 *
 * The RBF is rotation-invariant and minimizes bending energy, making it
 * ideal for smooth spatial interpolation in image differencing.
 *
 * @param r Euclidean distance between two points.
 * @return RBF kernel value φ(r).
 *
 * @note Special case: φ(0) = 0 by convention.
 */
static inline double tps_kernel(double r) {
  if (r < ZEROVAL) return 0.0;
  return r * r * log(r);
}

/**
 * @brief Assemble the augmented RBF matrix for TPS fitting.
 *
 * @details Solves the linear system:
 *   [ Φ  P ] [ w ]   [ c ]
 *   [ P' 0 ] [ v ] = [ 0 ]
 *
 * where:
 *   Φ[i,j] = φ(||pos[i] - pos[j]||)  (RBF kernel matrix)
 *   P[i,0:d] = polynomial basis (1, x, y) for linear trend
 *   w = RBF weights (N values per kernel component)
 *   v = polynomial coefficients (3 values: const, cx, cy)
 *   c = data values at control points
 *
 * The system is augmented with polynomial terms to ensure the interpolant
 * has a well-defined trend and avoids ill-conditioning.
 *
 * @param[in]  positions      (N, 2) array of (x, y) control point positions
 * @param[in]  n_points       Number of control points (N)
 * @param[out] matrix         (N+3, N+3) augmented RBF matrix (row-major)
 *
 * @note The caller must allocate matrix of size (n_points+3)²
 */
static void tps_assemble_matrix(double* positions, int n_points,
                                double* matrix) {
  int i, j, mat_size;
  double dx, dy, r;

  mat_size = n_points + 3;

  /* Fill RBF kernel block Φ (top-left) */
  for (i = 0; i < n_points; i++) {
    for (j = 0; j < n_points; j++) {
      dx = positions[2 * i] - positions[2 * j];
      dy = positions[2 * i + 1] - positions[2 * j + 1];
      r = sqrt(dx * dx + dy * dy);
      matrix[i * mat_size + j] = tps_kernel(r);
    }
  }

  /* Fill polynomial basis block P (top-right, n_points × 3) */
  for (i = 0; i < n_points; i++) {
    matrix[i * mat_size + n_points] = 1.0;                  /* const term */
    matrix[i * mat_size + n_points + 1] = positions[2 * i]; /* x term */
    matrix[i * mat_size + n_points + 2] = positions[2 * i + 1]; /* y term */
  }

  /* Fill transpose of polynomial block P' (bottom-left, 3 × n_points) */
  for (j = 0; j < n_points; j++) {
    matrix[n_points * mat_size + j] = 1.0;
    matrix[(n_points + 1) * mat_size + j] = positions[2 * j];
    matrix[(n_points + 2) * mat_size + j] = positions[2 * j + 1];
  }

  /* Fill zero block (bottom-right, 3 × 3) */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      matrix[(n_points + i) * mat_size + (n_points + j)] = 0.0;
    }
  }
}

/**
 * @brief Fit TPS RBF surface to scattered control points for one kernel
 *        component.
 *
 * @details Solves the augmented linear system via LU decomposition (LAPACK).
 *
 * @param[in]  positions         (N, 2) array of stamp center positions (pixels)
 * @param[in]  n_points          Number of stamps (control points)
 * @param[in]  coefficient_values (N,) array of fitted coefficients c_i for one
 *                                kernel component
 * @param[out] weights            (N,) RBF weights for this component
 * @param[out] poly_coeffs        (3,) polynomial trend coefficients [const, cx,
 *                                cy]
 * @return 0 on success, non-zero on singular matrix or other error
 *
 * @note Allocates and frees internal working arrays; caller provides output
 *       arrays.
 */
static int tps_fit_coefficients(double* positions, int n_points,
                                double* coefficient_values, double* weights,
                                double* poly_coeffs) {
  int mat_size, i, info;
  double* matrix;
  double* rhs;
  int* pivots;

  mat_size = n_points + 3;

  /* Allocate working arrays */
  matrix = (double*)xcalloc(mat_size * mat_size, sizeof(double));
  rhs = (double*)xcalloc(mat_size, sizeof(double));
  pivots = (int*)xcalloc(mat_size, sizeof(int));

  /* Assemble augmented RBF matrix */
  tps_assemble_matrix(positions, n_points, matrix);

  /* Copy coefficient values into RHS vector; pad with zeros for polynomial
   * terms */
  for (i = 0; i < n_points; i++) {
    rhs[i] = coefficient_values[i];
  }
  for (i = n_points; i < mat_size; i++) {
    rhs[i] = 0.0;
  }

  /* Solve via LU factorization (LAPACK dgesv) */
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, mat_size, 1, matrix, mat_size, pivots,
                       rhs, 1);

  if (info != 0) {
    LOG_ERROR("TPS fit failed: singular matrix (info=%d)", info);
    free(matrix);
    free(rhs);
    free(pivots);
    return 1;
  }

  /* Extract solution: RBF weights and polynomial coefficients */
  for (i = 0; i < n_points; i++) {
    weights[i] = rhs[i];
  }
  for (i = 0; i < 3; i++) {
    poly_coeffs[i] = rhs[n_points + i];
  }

  free(matrix);
  free(rhs);
  free(pivots);

  return 0;
}

/**
 * @brief Evaluate fitted TPS surface at an arbitrary position.
 *
 * @details Computes:
 *   c_i(x,y) = Σ_j w_j · φ(||(x,y) - pos_j||) + p[0] + p[1]·x + p[2]·y
 *
 * where φ is the TPS RBF kernel and p[] are polynomial trend coefficients.
 *
 * @param[in] eval_x, eval_y    Position at which to evaluate
 * @param[in] positions          (N, 2) array of control point positions
 * @param[in] n_points           Number of control points
 * @param[in] weights            (N,) RBF weights
 * @param[in] poly_coeffs        (3,) polynomial coefficients [const, cx, cy]
 * @return Interpolated coefficient value at (eval_x, eval_y)
 */
static double tps_evaluate(double eval_x, double eval_y, double* positions,
                           int n_points, double* weights,
                           double* poly_coeffs) {
  int i;
  double value, dx, dy, r;

  /* Initialize with polynomial trend */
  value = poly_coeffs[0] + poly_coeffs[1] * eval_x + poly_coeffs[2] * eval_y;

  /* Add RBF contributions */
  for (i = 0; i < n_points; i++) {
    dx = eval_x - positions[2 * i];
    dy = eval_y - positions[2 * i + 1];
    r = sqrt(dx * dx + dy * dy);
    value += weights[i] * tps_kernel(r);
  }

  return value;
}

/* =====================================================================
   TPS kernelSol Layout Management
   ===================================================================== */

/**
 * @brief Compute the size of kernelSol array needed for polynomial mode.
 *
 * Layout:
 *   [0]: reserved
 *   [1..nComp]: kernel polynomial coefficients
 *   [nComp+1..nComp+nbg_vec]: background polynomial coefficients
 *
 * where nComp = (nCompKer-1) * (kerOrder+1)*(kerOrder+2)/2
 *       nbg_vec = (bgOrder+1)*(bgOrder+2)/2
 */
static int kernelSol_size_polynomial(void) {
  int ncomp1 = nCompKer - 1;
  int ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  int ncomp = ncomp1 * ncomp2;
  int nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
  return ncomp + nbg_vec + 1;
}

/**

/**
 * @brief Get offset in kernelSol for RBF weights of kernel component.
 *
 * @param comp_idx Kernel component index (0..nCompKer-1)
 * @param n_stamps Number of stamps
 * @return Offset in kernelSol; access weights via kernelSol[offset..offset+n_stamps)
 */
static int kernelSol_offset_tps_weights(int comp_idx, int n_stamps) {
  int poly_size = kernelSol_size_polynomial();
  return poly_size + comp_idx * n_stamps;
}

/**
 * @brief Get offset in kernelSol for polynomial trend of kernel component.
 *
 * @param comp_idx Kernel component index (0..nCompKer-1)
 * @param n_stamps Number of stamps
 * @return Offset in kernelSol; access 3 poly coeffs via kernelSol[offset..offset+3)
 */
static int kernelSol_offset_tps_poly(int comp_idx, int n_stamps) {
  int poly_size = kernelSol_size_polynomial();
  return poly_size + nCompKer * n_stamps + comp_idx * 3;
}

/**
 * @brief Get offset in kernelSol for stamp positions.
 *
 * @param n_stamps Number of stamps
 * @return Offset in kernelSol; positions are stored as x0,y0,x1,y1,...
 */
static int kernelSol_offset_tps_positions(int n_stamps) {
  int poly_size = kernelSol_size_polynomial();
  return poly_size + nCompKer * n_stamps + 3 * nCompKer;
}

/**
 * @brief Get offset in kernelSol for RBF weights of background.
 *
 * @param n_stamps Number of stamps
 * @return Offset in kernelSol; access weights via kernelSol[offset..offset+n_stamps)
 */
static int kernelSol_offset_tps_bg_weights(int n_stamps) {
  int poly_size = kernelSol_size_polynomial();
  return poly_size + nCompKer * n_stamps + 3 * nCompKer + 2 * n_stamps;
}

/**
 * @brief Get offset in kernelSol for polynomial trend of background.
 *
 * @param n_stamps Number of stamps
 * @return Offset in kernelSol; access 3 poly coeffs via kernelSol[offset..offset+3)
 */
static int kernelSol_offset_tps_bg_poly(int n_stamps) {
  int poly_size = kernelSol_size_polynomial();
  return poly_size + nCompKer * n_stamps + 3 * nCompKer + 2 * n_stamps + n_stamps;
}

/**
 * @brief Fit thin plate spline surfaces to kernel coefficients after polynomial
 * solve.
 *
 * @details Called from fitKernel() after LAPACK solve (if useTPS==1). For each
 * kernel component:
 *   1. Evaluate the fitted polynomial at each stamp center
 *   2. Call tps_fit_coefficients() to fit RBF surface
 *   3. Store RBF weights and poly trend in kernelSol
 *
 * @param[in]  stamps     Array of stamps with centers
 * @param[in]  n_stamps   Total number of stamps
 * @param[in,out] kernelSol  Solution vector; polynomial part is read, TPS part
 *                           is filled
 * @return 0 on success, non-zero if any TPS fit fails
 */
static int tps_fit_kernel(stamp_struct* stamps, int n_stamps,
                          double* kernelSol) {
  int comp_idx, stamp_idx, i;
  int ncomp1, ncomp2, ncomp, solutionIdx;
  double *stamp_positions, *poly_coeffs_at_stamps;
  double *tps_weights, *tps_poly;
  double normalizedX, normalizedY, polyBasisX, polyBasisY;
  double halfPixX, halfPixY;
  int polyDegX, polyDegY;

  ncomp1 = nCompKer - 1;
  ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  ncomp = ncomp1 * ncomp2;

  halfPixX = 0.5 * rPixX;
  halfPixY = 0.5 * rPixY;

  LOG_DEBUG("TPS refitting: %d stamps, %d kernel components", n_stamps, nCompKer);

  /* Allocate working arrays */
  stamp_positions = (double*)xcalloc(2 * n_stamps, sizeof(double));
  poly_coeffs_at_stamps = (double*)xcalloc(n_stamps, sizeof(double));

  /* Extract stamp centers (in image pixel coordinates) */
  for (stamp_idx = 0; stamp_idx < n_stamps; stamp_idx++) {
    stamp_positions[2 * stamp_idx] = (double)stamps[stamp_idx].x;
    stamp_positions[2 * stamp_idx + 1] = (double)stamps[stamp_idx].y;
  }

  /* Store stamp positions in kernelSol for later evaluation */
  {
    int pos_offset = kernelSol_offset_tps_positions(n_stamps);
    for (i = 0; i < 2 * n_stamps; i++) {
      kernelSol[pos_offset + i] = stamp_positions[i];
    }
  }

  /* For each kernel component (excluding constant term), fit TPS surface */
  for (comp_idx = 1; comp_idx < nCompKer; comp_idx++) {
    solutionIdx = 2 + (comp_idx - 1) * ncomp2;

    /* Evaluate polynomial at each stamp to get per-stamp coefficients */
    for (stamp_idx = 0; stamp_idx < n_stamps; stamp_idx++) {
      normalizedX =
          (stamp_positions[2 * stamp_idx] - halfPixX) / halfPixX;
      normalizedY =
          (stamp_positions[2 * stamp_idx + 1] - halfPixY) / halfPixY;

      poly_coeffs_at_stamps[stamp_idx] = 0.0;
      polyBasisX = 1.0;
      for (polyDegX = 0; polyDegX <= kerOrder; polyDegX++) {
        polyBasisY = 1.0;
        for (polyDegY = 0; polyDegY <= kerOrder - polyDegX; polyDegY++) {
          poly_coeffs_at_stamps[stamp_idx] +=
              kernelSol[solutionIdx++] * polyBasisX * polyBasisY;
          polyBasisY *= normalizedY;
        }
        polyBasisX *= normalizedX;
      }
    }

    /* Fit TPS surface to per-stamp coefficients */
    {
      int weights_offset = kernelSol_offset_tps_weights(comp_idx, n_stamps);
      int poly_offset = kernelSol_offset_tps_poly(comp_idx, n_stamps);
      int ret;

      tps_weights = kernelSol + weights_offset;
      tps_poly = kernelSol + poly_offset;

      ret = tps_fit_coefficients(stamp_positions, n_stamps,
                                 poly_coeffs_at_stamps, tps_weights, tps_poly);
      if (ret != 0) {
        LOG_ERROR("TPS fit failed for kernel component %d", comp_idx);
        free(stamp_positions);
        free(poly_coeffs_at_stamps);
        return 1;
      }

      LOG_DEBUG("TPS component %d: fitted %d RBF weights, 3 poly trends", comp_idx,
                n_stamps);
    }
  }

  free(stamp_positions);
  free(poly_coeffs_at_stamps);

  return 0;
}

/**
 * @brief Fit thin plate spline surface to background coefficients after polynomial
 * solve.
 *
 * @details Called from fitKernel() after tps_fit_kernel() (if useTPS==1). Evaluates
 * the fitted background polynomial at each stamp center and fits an RBF surface.
 *
 * @param[in]  stamps     Array of stamps with centers
 * @param[in]  n_stamps   Total number of stamps
 * @param[in,out] kernelSol  Solution vector; background polynomial part is read,
 *                           background TPS part is filled
 * @return 0 on success, non-zero if TPS fit fails
 */
static int tps_fit_background(stamp_struct* stamps, int n_stamps,
                              double* kernelSol) {
  int stamp_idx, i;
  double *stamp_positions, *bg_coeffs_at_stamps;
  double *tps_weights, *tps_poly;
  double normalizedX, normalizedY, polyBasisX, polyBasisY;
  double halfPixX, halfPixY;
  int polyDegX, polyDegY, bgSolutionIdx;
  int nbg_vec;

  nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;

  halfPixX = 0.5 * rPixX;
  halfPixY = 0.5 * rPixY;

  LOG_DEBUG("TPS background refitting: %d stamps", n_stamps);

  /* Allocate working arrays */
  stamp_positions = (double*)xcalloc(2 * n_stamps, sizeof(double));
  bg_coeffs_at_stamps = (double*)xcalloc(n_stamps, sizeof(double));

  /* Extract stamp centers (in image pixel coordinates) */
  for (stamp_idx = 0; stamp_idx < n_stamps; stamp_idx++) {
    stamp_positions[2 * stamp_idx] = (double)stamps[stamp_idx].x;
    stamp_positions[2 * stamp_idx + 1] = (double)stamps[stamp_idx].y;
  }

  /* Compute background polynomial offset in kernelSol */
  {
    int ncomp1 = nCompKer - 1;
    int ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    int ncomp = ncomp1 * ncomp2;
    bgSolutionIdx = ncomp + 1;
  }

  /* Evaluate background polynomial at each stamp */
  for (stamp_idx = 0; stamp_idx < n_stamps; stamp_idx++) {
    normalizedX = (stamp_positions[2 * stamp_idx] - halfPixX) / halfPixX;
    normalizedY = (stamp_positions[2 * stamp_idx + 1] - halfPixY) / halfPixY;

    bg_coeffs_at_stamps[stamp_idx] = 0.0;
    polyBasisX = 1.0;
    for (polyDegX = 0; polyDegX <= bgOrder; polyDegX++) {
      polyBasisY = 1.0;
      for (polyDegY = 0; polyDegY <= bgOrder - polyDegX; polyDegY++) {
        bg_coeffs_at_stamps[stamp_idx] +=
            kernelSol[bgSolutionIdx++] * polyBasisX * polyBasisY;
        polyBasisY *= normalizedY;
      }
      polyBasisX *= normalizedX;
    }
  }

  /* Fit TPS surface to per-stamp background coefficients */
  {
    int weights_offset = kernelSol_offset_tps_bg_weights(n_stamps);
    int poly_offset = kernelSol_offset_tps_bg_poly(n_stamps);
    int ret;

    tps_weights = kernelSol + weights_offset;
    tps_poly = kernelSol + poly_offset;

    ret = tps_fit_coefficients(stamp_positions, n_stamps,
                               bg_coeffs_at_stamps, tps_weights, tps_poly);
    if (ret != 0) {
      LOG_ERROR("TPS fit failed for background");
      free(stamp_positions);
      free(bg_coeffs_at_stamps);
      return 1;
    }

    LOG_DEBUG("TPS background: fitted %d RBF weights, 3 poly trends", n_stamps);
  }

  free(stamp_positions);
  free(bg_coeffs_at_stamps);

  return 0;
}

/* =====================================================================
 */

/**
 * @brief Precompute all kernel basis function images and store them in
 * kernel_vec.
 *
 * @details Iterates over every Gaussian component and its associated polynomial
 * orders (deg_fixe[ig]), calling kernel_vector() for each (ig, idegx, idegy)
 * triplet.  The resulting fwKernel×fwKernel images are the fixed basis profiles
 * φᵢ in the kernel expansion K(x,y) = Σᵢ cᵢ(x,y)·φᵢ (Alard & Lupton 1998,
 * Eq. 2).  This function is called exactly once during initialisation; the
 * populated global array kernel_vec is subsequently used by every convolution
 * and model-evaluation call.
 */
void getKernelVec() {
  /* Dispatch basis initialization based on iBasisType */
  if (iBasisType == BASIS_TYPE_DELTA) {
    LOG_ERROR("Delta function basis not yet implemented (Phase 2+)");
    /* TODO: int nBasisFuncs = 0;
           init_delta_basis_grid(&nBasisFuncs);
           nCompKer = nBasisFuncs; */
    return;
  }

  /* Gaussian basis initialization (default) */
  int gaussIdx, idegx, idegy, nvec;
  int ren;

  nvec = 0;
  for (gaussIdx = 0; gaussIdx < ngauss; gaussIdx++) {
    for (idegx = 0; idegx <= deg_fixe[gaussIdx]; idegx++) {
      for (idegy = 0; idegy <= deg_fixe[gaussIdx] - idegx; idegy++) {
        /* stores kernel weight mask for each order */
        kernel_vec[nvec] = kernel_vector(nvec, idegx, idegy, gaussIdx, &ren);
        nvec++;
      }
    }
  }
}

/**
 * @brief Populate a stamp's convolved-image vectors and build its per-stamp
 *        normal-equation matrices.
 *
 * @details For the current substamp (stamp->sscnt):
 *  1. Calls xy_conv_stamp() for each kernel basis function to fill
 *     stamp->vectors[0..nCompKer-1] with the image convolved by that basis
 *     profile.
 *  2. Fills stamp->vectors[nCompKer..] with polynomial background basis
 *     images x^i · y^j evaluated over the substamp pixel grid (normalised to
 *     [-1, 1]).
 *  3. Calls build_matrix0() and build_scprod0() to accumulate the per-stamp
 *     least-squares cross-product matrices that will later be combined by
 *     build_matrix() and build_scprod().
 *
 * If all substamps for this stamp have been exhausted (sscnt >= nss) the
 * stamp is rejected and the function returns 1.
 *
 * @param stamp   Pointer to the stamp being processed; stamp->sscnt selects
 *                the active substamp.
 * @param imConv  Flat pixel array of the image to be convolved (the "template"
 *                or "science" image depending on the convolution direction).
 * @param imRef   Flat pixel array of the reference image (the target of the
 *                kernel fit, used to build scprod).
 * @return 0 on success, 1 if the stamp is rejected (out of substamps or bad
 *         substamp data).
 */
int fillStamp(stamp_struct* stamp, float* imConv, float* imRef) {
  int renormalizeFlag = 0;
  int pixelX, pixelY, stampCenterX, stampCenterY, pixelOffsetX, pixelOffsetY,
      idegx, idegy, bgDegX, bgDegY, substampIndexX, substampIndexY,
      vectorCompIdx, gaussianCompIdx, vectorComponentIdx;
  double polyBasisX, polyBasisY, normalizedX, normalizedY;
  double* im;
  float halfPixX, halfPixY;

  halfPixX = 0.5 * rPixX;
  halfPixY = 0.5 * rPixY;

  LOG_DEBUG_INDENT(4, "xs  : %4i ys  : %4i sig: %6.3f sscnt: %4i nss: %4i",
                   stamp->x, stamp->y, stamp->chi2, stamp->sscnt, stamp->nss);
  if (stamp->sscnt >= stamp->nss) {
    /* have gone through all the good substamps, reject this stamp */
    /*if (verbose >= 2) fprintf(stderr, "    ******** REJECT stamp (out of
     * substamps)\n");*/
    LOG_DEBUG_INDENT(4, "Stamp rejected");
    return 1;
  }

  /* ====================================================================
     POLYNOMIAL BASIS CONSTRUCTION
     ====================================================================
     Iterate over kernel basis triplets (gaussianCompIdx, idegx, idegy):
     - gaussianCompIdx: Gaussian component index [0, ngauss)
     - idegx, idegy: polynomial degrees such that idegx + idegy <=
     deg_fixe[gaussianCompIdx]

     This implements the triangular basis expansion:
       K(x,y) = Sum_i c_i(x,y) * phi_i
     where phi_i = exp(-(x^2+y^2)*sigma^2) * x^deg_x * y^deg_y

     The counter vectorComponentIdx sequentially assigns indices to each
     (gaussianCompIdx, idegx, idegy) triplet; these indices are used by
     xy_conv_stamp() to access precomputed convolution responses stored in
     stamp->vectors[vectorComponentIdx].

     Alard & Lupton (1998), Sect. 2: "The kernel is expanded as a sum of
     Gaussian basis elements with polynomial spatial weighting."
     Reference: https://iopscience.iop.org/article/10.1086/305984
     ==================================================================== */
  vectorComponentIdx = 0;
  for (gaussianCompIdx = 0; gaussianCompIdx < ngauss; gaussianCompIdx++) {
    for (idegx = 0; idegx <= deg_fixe[gaussianCompIdx]; idegx++) {
      for (idegy = 0; idegy <= deg_fixe[gaussianCompIdx] - idegx; idegy++) {
        renormalizeFlag = 0;
        /* Detect odd-degree terms: (deg/2)*2 - deg evaluates to
           -1 if deg is odd, 0 if even. Used to trigger renormalization
           of even-parity basis functions (those with pixelOffsetX==0 &&
           pixelOffsetY==0). */
        pixelOffsetX = (idegx / 2) * 2 - idegx;
        pixelOffsetY = (idegy / 2) * 2 - idegy;
        if (pixelOffsetX == 0 && pixelOffsetY == 0 && vectorComponentIdx > 0)
          renormalizeFlag = 1; /* Set renormalization flag for higher-order
                                  even-parity basis */

        /* fill stamp->vectors[vectorComponentIdx] with convolved image */
        /* image is convolved with functional form of kernel, fit later for
         * amplitude */
        xy_conv_stamp(stamp, imConv, vectorComponentIdx, renormalizeFlag);
        ++vectorComponentIdx;
      }
    }
  }

  /* get the krefArea data */
  if (cutSStamp(stamp, imRef)) return 1;

  /* fill stamp->vectors[vectorComponentIdx+++] with x^(bg) * y^(bg) for
   * background fit */
  stampCenterX = stamp->xss[stamp->sscnt];
  stampCenterY = stamp->yss[stamp->sscnt];
  substampIndexX = stampCenterX - hwKSStamp;
  substampIndexY = stampCenterY - hwKSStamp;
  for (pixelX = stampCenterX - hwKSStamp; pixelX <= stampCenterX + hwKSStamp;
       pixelX++) {
    normalizedX = (pixelX - halfPixX) / halfPixX;

    for (pixelY = stampCenterY - hwKSStamp; pixelY <= stampCenterY + hwKSStamp;
         pixelY++) {
      /* fprintf(stderr, "%d %d %d %d %d %d\n", k, stampCenterX,
       * stampCenterY,pixelX, pixelY, fwKSStamp); */
      normalizedY = (pixelY - halfPixY) / halfPixY;

      polyBasisX = 1.0;
      vectorCompIdx = vectorComponentIdx;
      for (bgDegX = 0; bgDegX <= bgOrder; bgDegX++) {
        polyBasisY = 1.0;
        for (bgDegY = 0; bgDegY <= bgOrder - bgDegX; bgDegY++) {
          im = stamp->vectors[vectorCompIdx];
          im[pixelX - substampIndexX + fwKSStamp * (pixelY - substampIndexY)] =
              polyBasisX * polyBasisY;
          polyBasisY *= normalizedY;
          ++vectorCompIdx;
        }
        polyBasisX *= normalizedX;
      }
    }
  }

  /* build stamp->mat from stamp->vectors */
  build_matrix0(stamp);
  /* build stamp->scprod from stamp->vectors and imRef */
  build_scprod0(stamp, imRef);

  return 0;
}

/**
 * @brief Compute the fwKernel×fwKernel image for one Gaussian-polynomial kernel
 *        basis function.
 *
 * @details The n-th basis function is the outer product of two 1-D filters:
 *   filter_x[k] = exp(-x² · sigma_gauss[ig]) · x^deg_x
 *   filter_y[k] = exp(-x² · sigma_gauss[ig]) · x^deg_y
 * where x runs over [-hwKernel, +hwKernel].
 *
 * For even-parity bases (deg_x and deg_y both even) the filters are
 * normalised so that they sum to unity.  Additionally, for n > 0 the zeroth
 * basis image is subtracted, making higher-order terms represent differential
 * corrections to the zeroth-order (pure Gaussian) kernel as described in
 * Alard & Lupton (1998), §2.  When this subtraction is applied *ren is set to
 * 1 (the "renormalise" flag propagated to xy_conv_stamp).
 *
 * @param n     Sequential index of this basis function within kernel_vec.
 * @param deg_x Polynomial degree of the x-direction filter.
 * @param deg_y Polynomial degree of the y-direction filter.
 * @param ig    Gaussian component index (selects sigma_gauss[ig]).
 * @param ren   Output flag: set to 1 when the zeroth basis is subtracted.
 * @return Pointer to a newly allocated fwKernel*fwKernel double array; the
 *         caller stores this in kernel_vec[n] and must eventually free it.
 */
double* kernel_vector(int n, int deg_x, int deg_y, int ig, int* ren) {
  double *vector = NULL, *kernel0 = NULL;
  int kernelRow, kernelCol, xFilterIdx, xParityCheck, yParityCheck, xKernelIdx,
      pixelIdx;
  double filterXSum, filterYSum, kernelPosition, gaussianVal;

  vector = (double*)malloc(fwKernel * fwKernel * sizeof(double));
  /* Detect parity: (deg/2)*2 - deg = 0 if even, -1 if odd.
     Odd-degree filters (xParityCheck != 0 || yParityCheck != 0) retain full
     normalization; even-degree filters (xParityCheck==0 && yParityCheck==0) are
     normalized and (for n>0) have the n=0 basis subtracted to form differential
     corrections. */
  xParityCheck = (deg_x / 2) * 2 - deg_x;
  yParityCheck = (deg_y / 2) * 2 - deg_y;
  filterXSum = filterYSum = 0.0;
  *ren = 0;

  for (xKernelIdx = 0; xKernelIdx < fwKernel; xKernelIdx++) {
    kernelPosition = (double)(xKernelIdx - hwKernel);
    xFilterIdx = xKernelIdx + n * fwKernel;
    gaussianVal = exp(-kernelPosition * kernelPosition * sigma_gauss[ig]);
    filter_x[xFilterIdx] = gaussianVal * pow(kernelPosition, deg_x);
    filter_y[xFilterIdx] = gaussianVal * pow(kernelPosition, deg_y);
    filterXSum += filter_x[xFilterIdx];
    filterYSum += filter_y[xFilterIdx];
  }

  if (n > 0) kernel0 = kernel_vec[0];

  filterXSum = 1. / filterXSum;
  filterYSum = 1. / filterYSum;

  if (xParityCheck == 0 && yParityCheck == 0) {
    for (xKernelIdx = 0; xKernelIdx < fwKernel; xKernelIdx++) {
      filter_x[xKernelIdx + n * fwKernel] *= filterXSum;
      filter_y[xKernelIdx + n * fwKernel] *= filterYSum;
    }

    for (kernelRow = 0; kernelRow < fwKernel; kernelRow++) {
      for (kernelCol = 0; kernelCol < fwKernel; kernelCol++) {
        vector[kernelRow + fwKernel * kernelCol] =
            filter_x[kernelRow + n * fwKernel] *
            filter_y[kernelCol + n * fwKernel];
      }
    }

    if (n > 0) {
      for (pixelIdx = 0; pixelIdx < fwKernel * fwKernel; pixelIdx++) {
        vector[pixelIdx] -= kernel0[pixelIdx];
      }
      *ren = 1;
    }
  } else {
    for (kernelRow = 0; kernelRow < fwKernel; kernelRow++) {
      for (kernelCol = 0; kernelCol < fwKernel; kernelCol++) {
        vector[kernelRow + fwKernel * kernelCol] =
            filter_x[kernelRow + n * fwKernel] *
            filter_y[kernelCol + n * fwKernel];
      }
    }
  }
  return vector;
}

/**
 * @brief Convolve the substamp region of an image with the n-th separable
 *        Gaussian-polynomial kernel basis function.
 *
 * @details **Separability optimization:** Gaussian-polynomial filters are
 * separable: φ(x,y) = filter_x[x] * filter_y[y] where filter_x[k] =
 * exp(-k^2*σ^2) * k^deg_x and similarly for filter_y.
 *
 * Standard 2D convolution costs O(k^2*n^2) FLOPs (k = kernel width, n = image
 * dimension). Separable 1D passes reduce this to O(2*k*n^2), a speedup factor
 * of k/2 (typically 5–10×). This is the primary reason HOTPANTS achieves
 * real-time performance on large images.
 *
 * The implementation uses the separability of the Gaussian-polynomial filter to
 * perform the convolution in two 1-D passes (first along y, then along x),
 * storing the fwKSStamp×fwKSStamp result in stamp->vectors[n].  This is the
 * performance- critical inner loop (~60 % of total runtime according to
 * profiling).
 *
 * If the renormalise flag @p ren is set (i.e. kernel_vector() returned ren=1
 * for this basis), the zeroth-basis convolved image stamp->vectors[0] is
 * subtracted pixel-by-pixel from the result, matching the differential
 * construction of the higher-order basis images.
 *
 * Reference: Alard & Lupton (1998), Sect. 2; classic signal-processing
 * optimization.
 *
 * @param stamp  Stamp whose current substamp (stamp->sscnt) defines the pixel
 *               region to convolve.
 * @param image  Full-frame input image to be convolved.
 * @param n      Index of the kernel basis function; selects filter_x/filter_y
 *               rows and the output slot stamp->vectors[n].
 * @param ren    If non-zero, subtract stamp->vectors[0] from the result after
 *               convolution (renormalisation step for higher-order bases).
 */
void xy_conv_stamp(stamp_struct* stamp, float* image, int n, int ren) {
  int stampColIdx, stampRowIdx, kerOffsetX, kerOffsetY, xij, sub_width, xi, yi,
      imgColIdx, imgRowIdx, pixelIdx;
  double *v0, *imc;

  xi = stamp->xss[stamp->sscnt];
  yi = stamp->yss[stamp->sscnt];
  imc = stamp->vectors[n];

  sub_width = fwKSStamp + fwKernel - 1;

  /* pull area to convolve out of full reference image region */
  /* convolve with y filter */
  for (imgColIdx = xi - hwKSStamp - hwKernel;
       imgColIdx <= xi + hwKSStamp + hwKernel; imgColIdx++) {
    for (imgRowIdx = yi - hwKSStamp; imgRowIdx <= yi + hwKSStamp; imgRowIdx++) {
      xij = imgColIdx - xi + sub_width / 2 +
            sub_width * (imgRowIdx - yi + hwKSStamp);
      temp[xij] = 0.0;
      for (kerOffsetY = -hwKernel; kerOffsetY <= hwKernel; kerOffsetY++) {
        temp[xij] += image[imgColIdx + rPixX * (imgRowIdx + kerOffsetY)] *
                     filter_y[hwKernel - kerOffsetY + n * fwKernel];
      }
    }
  }

  /* convolve with x filter */
  for (stampRowIdx = -hwKSStamp; stampRowIdx <= hwKSStamp; stampRowIdx++) {
    for (stampColIdx = -hwKSStamp; stampColIdx <= hwKSStamp; stampColIdx++) {
      xij = stampColIdx + hwKSStamp + fwKSStamp * (stampRowIdx + hwKSStamp);
      imc[xij] = 0.0;
      for (kerOffsetX = -hwKernel; kerOffsetX <= hwKernel; kerOffsetX++) {
        imc[xij] += temp[stampColIdx + kerOffsetX + sub_width / 2 +
                         sub_width * (stampRowIdx + hwKSStamp)] *
                    filter_x[hwKernel - kerOffsetX + n * fwKernel];
      }
    }
  }

  if (ren) {
    v0 = stamp->vectors[0];
    for (pixelIdx = 0; pixelIdx < fwKSStamp * fwKSStamp; pixelIdx++)
      imc[pixelIdx] -= v0[pixelIdx];
  }

  return;
}


/**
 * @brief Solve A·x = b in-place where A is symmetric positive-definite.
 *
 * @details Adaptor between the Numerical Recipes 1-based double** convention
 * used throughout this file and the LAPACKE C interface.  Copies the matrix and
 * RHS into contiguous 0-based row-major buffers, runs Cholesky factorisation
 * (dpotrf) followed by the triangular solve (dpotrs), then writes the solution
 * back into b[1..n].  Cholesky is the correct choice here because the
 * normal-equations matrix M = Σᵢ AᵢᵀAᵢ is symmetric positive-definite by
 * construction; it is both faster and more numerically stable than LU.
 *
 * @param a  1-indexed n×n SPD matrix a[1..n][1..n]; modified in place during
 *           factorisation (upper triangle overwritten with Cholesky factor).
 * @param n  Matrix dimension.
 * @param b  1-indexed RHS vector b[1..n]; overwritten with the solution x on
 *           successful return.
 * @return 0 on success, 1 if the matrix is not positive-definite or allocation
 *         fails.
 */
static int solve_spd(double** a, int n, double* b) {
  lapack_int solveStatus;
  int matrixRowIdx, matrixColIdx;
  double* flatMatrix = (double*)malloc((size_t)n * n * sizeof(double));
  double* rhsVector = (double*)malloc((size_t)n * sizeof(double));
  if (!flatMatrix || !rhsVector) {
    free(flatMatrix);
    free(rhsVector);
    return 1;
  }

  for (matrixRowIdx = 0; matrixRowIdx < n; matrixRowIdx++) {
    for (matrixColIdx = 0; matrixColIdx < n; matrixColIdx++)
      flatMatrix[matrixRowIdx * n + matrixColIdx] =
          a[matrixRowIdx + 1][matrixColIdx + 1];
    rhsVector[matrixRowIdx] = b[matrixRowIdx + 1];
  }

  solveStatus = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', (lapack_int)n, flatMatrix,
                               (lapack_int)n);
  if (solveStatus == 0)
    solveStatus = LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'U', (lapack_int)n, 1,
                                 flatMatrix, (lapack_int)n, rhsVector, 1);
  if (solveStatus == 0)
    for (matrixRowIdx = 0; matrixRowIdx < n; matrixRowIdx++)
      b[matrixRowIdx + 1] = rhsVector[matrixRowIdx];
  else
    LOG_WARNING(
        "solve_spd failed (LAPACKE status=%d); "
        "kernel solution may be unreliable",
        (int)solveStatus);

  free(flatMatrix);
  free(rhsVector);
  return (solveStatus != 0) ? 1 : 0;
}

/**
 * @brief Solve for the spatially-varying kernel by fitting the global normal
 *        equations, iterating until stamp sigma-clipping converges.
 *
 * @details Implements the full Alard & Lupton (1998) kernel solution:
 *  1. Assembles the global normal-equations matrix M via build_matrix() and
 *     the right-hand-side vector b via build_scprod().
 *  2. Solves M·a = b via LAPACK Cholesky decomposition (dpotrf + dpotrs);
 *     the solution vector a is written into @p kernelSol.
 *  3. Calls check_again() to evaluate each stamp's residual sigma and reject
 *     or replace bad substamps; if any stamp is flagged the matrices are
 *     rebuilt and the linear system is solved again.
 *  4. Iterates steps 1–3 until check_again() reports no further changes.
 *
 * The normal-equations matrix is symmetric positive-definite; Cholesky
 * decomposition is optimal for this structure.
 *
 * @param stamps              Array of nS stamps used for the fit.
 * @param imRef               Full-frame reference image (the fit target).
 * @param imConv              Full-frame image to be convolved (the template or
 *                            science image).
 * @param imNoise             Per-pixel noise (sigma) image used when evaluating
 *                            stamp residuals inside check_again().
 * @param kernelSol           Output: solution vector of length nCompTotal+1
 *                            containing the spatially-varying kernel and
 *                            background coefficients.
 * @param meansigSubstamps    Output: sigma-clipped mean residual over stamps
 *                            at convergence.
 * @param scatterSubstamps    Output: sigma-clipped scatter of residuals over
 *                            stamps at convergence.
 * @param NskippedSubstamps   Output: number of substamps excluded from the
 *                            final solution.
 */
void fitKernel(stamp_struct* stamps, float* imRef, float* imConv,
               float* imNoise, double* kernelSol, double* meansigSubstamps,
               double* scatterSubstamps, int* NskippedSubstamps) {
  double** matrix;
  char convergedFlag;
  int matrixRowIdx, mat_size;
  int ncomp1, ncomp2, ncomp, nbg_vec;

  ncomp1 = nCompKer - 1;
  ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  ncomp = ncomp1 * ncomp2;
  nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;

  mat_size = ncomp1 * ncomp2 + nbg_vec + 1;

  LOG_DEBUG("Mat_size: %i ncomp2: %i ncomp1: %i nbg_vec: %i", mat_size, ncomp2,
            ncomp1, nbg_vec);

  /* allocate fitting matrix */
  matrix = (double**)malloc((mat_size + 1) * sizeof(double*));
  for (matrixRowIdx = 0; matrixRowIdx <= mat_size; matrixRowIdx++)
    matrix[matrixRowIdx] = (double*)malloc((mat_size + 1) * sizeof(double));

  /* allocate weight matrix */
  wxy = (double**)malloc(nS * sizeof(double*));
  for (matrixRowIdx = 0; matrixRowIdx < nS; matrixRowIdx++)
    wxy[matrixRowIdx] = (double*)malloc(ncomp2 * sizeof(double));

  LOG_DEBUG("Expanding Matrix For Full Fit");
  build_matrix(stamps, nS, matrix);
  build_scprod(stamps, nS, imRef, kernelSol);

  solve_spd(matrix, mat_size, kernelSol);

  LOG_DEBUG("Checking again");
  convergedFlag =
      check_again(stamps, kernelSol, imConv, imRef, imNoise, meansigSubstamps,
                  scatterSubstamps, NskippedSubstamps);

  while (convergedFlag) {
    LOG_PROGRESS("Re-Expanding Matrix");
    build_matrix(stamps, nS, matrix);
    build_scprod(stamps, nS, imRef, kernelSol);

    solve_spd(matrix, mat_size, kernelSol);

    LOG_PROGRESS("Checking again");
    convergedFlag =
        check_again(stamps, kernelSol, imConv, imRef, imNoise, meansigSubstamps,
                    scatterSubstamps, NskippedSubstamps);
  }
  LOG_PROGRESS("Sigma clipping of bad stamps converged, kernel determined");

  /* If TPS mode is enabled, refit kernel and background coefficients using RBF interpolation */
  if (useTPS) {
    int ret = tps_fit_kernel(stamps, nS, kernelSol);
    if (ret != 0) {
      LOG_ERROR("TPS kernel refitting failed; using polynomial solution only");
      useTPS = 0; /* Fall back to polynomial evaluation */
    } else {
      LOG_PROGRESS("TPS kernel refitting complete; RBF surfaces fitted for %d components",
                   nCompKer);

      /* Also fit TPS surface to background */
      ret = tps_fit_background(stamps, nS, kernelSol);
      if (ret != 0) {
        LOG_ERROR("TPS background refitting failed; using polynomial background");
        /* Note: we don't disable useTPS here; kernel TPS still works, just background reverts to polynomial */
      } else {
        LOG_PROGRESS("TPS background refitting complete");
      }
    }
  }

  for (matrixRowIdx = 0; matrixRowIdx <= mat_size; matrixRowIdx++)
    free(matrix[matrixRowIdx]);
  for (matrixRowIdx = 0; matrixRowIdx < nS; matrixRowIdx++)
    free(wxy[matrixRowIdx]);
  free(matrix);
  free(wxy);

  return;
}

/**
 * @brief Accumulate the per-stamp least-squares cross-product matrix from the
 *        precomputed convolved-image vectors.
 *
 * @details Fills stamp->mat with the pixel-space inner products of the kernel
 * basis convolution images, corresponding to the Q matrix of Alard & Lupton
 * (1998), Eq. 3:
 *   stamp->mat[i][j] = Σ_k  vectors[i][k] · vectors[j][k]
 * for i,j in [1..nCompKer] (kernel basis terms) and the background coupling
 * terms (vectors[nCompKer]).  Only the upper triangle and the background row
 * are filled here; build_matrix() copies to the lower triangle when assembling
 * the global matrix.
 *
 * @param stamp  Stamp whose stamp->vectors have been filled by fillStamp();
 *               the result is stored in stamp->mat.
 */
void build_matrix0(stamp_struct* stamp) {
  int kernelBasisIdx1, kernelBasisIdx2, pixStamp, pixelIdx, kernelBasisIdx,
      bgVectorIdx = 0;
  int ncomp1, ncomp2, ncomp, nbg_vec;
  double p0, q;
  double** vec;

  ncomp1 = nCompKer;
  ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  ncomp = ncomp1 * ncomp2;
  nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;

  pixStamp = fwKSStamp * fwKSStamp;

  vec = stamp->vectors;

  /* loop over the convolved images created by xy_conv_stamp() */
  /* each level represents ngauss and deg_gauss */
  for (kernelBasisIdx1 = 0; kernelBasisIdx1 < ncomp1; kernelBasisIdx1++) {
    for (kernelBasisIdx2 = 0; kernelBasisIdx2 <= kernelBasisIdx1;
         kernelBasisIdx2++) {
      q = 0.0;
      /* integrate W_m1 and W_m2 (sum over all pixels) */
      for (pixelIdx = 0; pixelIdx < pixStamp; pixelIdx++)
        q += vec[kernelBasisIdx1][pixelIdx] * vec[kernelBasisIdx2][pixelIdx];

      /* Q from Eqn 3. in Alard */
      stamp->mat[kernelBasisIdx1 + 1][kernelBasisIdx2 + 1] = q;
    }
  }

  for (kernelBasisIdx = 0; kernelBasisIdx < ncomp1; kernelBasisIdx++) {
    /* bgVectorIdx = index into the background vector array (the constant
     * background term) */
    bgVectorIdx = ncomp1;

    p0 = 0.0;
    /* integrate convolved images and first order background (equals 1
     * everywhere!)*/
    for (pixelIdx = 0; pixelIdx < pixStamp; pixelIdx++)
      p0 += vec[kernelBasisIdx][pixelIdx] * vec[bgVectorIdx][pixelIdx];
    stamp->mat[ncomp1 + 1][kernelBasisIdx + 1] = p0;
  }

  /* integrate first order background with itself */
  /* NOTE : DON'T MASK K HERE - BACKGROUND! */
  for (pixelIdx = 0, q = 0.0; pixelIdx < pixStamp; pixelIdx++)
    q += vec[bgVectorIdx][pixelIdx] * vec[ncomp1][pixelIdx];
  stamp->mat[ncomp1 + 1][ncomp1 + 1] = q;

  return;
}

/**
 * @brief Compute the per-stamp right-hand-side vector (scalar products of
 *        basis convolutions with the reference image).
 *
 * @details Fills stamp->scprod with the pixel-space inner products
 *   stamp->scprod[i] = Σ_k  vectors[i][k] · image[xc+rPixX*(yc+yi)]
 * corresponding to the b vector in the normal equations M·a = b, implementing
 * Alard & Lupton (1998), Eq. 4.  The first nCompKer entries correspond to
 * kernel basis terms; the final entry is the cross product with the constant
 * background basis.
 *
 * @param stamp  Stamp whose stamp->vectors have been filled by fillStamp().
 * @param image  Full-frame reference image (fitting target).
 */
void build_scprod0(stamp_struct* stamp, float* image) {
  int stampColIdx, stampRowIdx, xi, yi, kernelBasisIdx, pixelIdx;
  int ncomp1, ncomp2, ncomp, nbg_vec;
  double p0, q;
  double** vec;

  ncomp1 = nCompKer;
  ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  ncomp = ncomp1 * ncomp2;
  nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;

  vec = stamp->vectors;
  xi = stamp->xss[stamp->sscnt];
  yi = stamp->yss[stamp->sscnt];

  /* Do eqn 4. in Alard */

  /* Multiply each order's convolved image with reference image */
  for (kernelBasisIdx = 0; kernelBasisIdx < ncomp1; kernelBasisIdx++) {
    p0 = 0.0;
    for (stampColIdx = -hwKSStamp; stampColIdx <= hwKSStamp; stampColIdx++) {
      for (stampRowIdx = -hwKSStamp; stampRowIdx <= hwKSStamp; stampRowIdx++) {
        pixelIdx =
            stampColIdx + hwKSStamp + fwKSStamp * (stampRowIdx + hwKSStamp);
        p0 += vec[kernelBasisIdx][pixelIdx] *
              image[stampColIdx + xi + rPixX * (stampRowIdx + yi)];
      }
    }
    stamp->scprod[kernelBasisIdx + 1] = p0;
  }

  /* Multiply first order background model with reference image */
  q = 0.0;
  for (stampColIdx = -hwKSStamp; stampColIdx <= hwKSStamp; stampColIdx++) {
    for (stampRowIdx = -hwKSStamp; stampRowIdx <= hwKSStamp; stampRowIdx++) {
      pixelIdx =
          stampColIdx + hwKSStamp + fwKSStamp * (stampRowIdx + hwKSStamp);
      q += vec[ncomp1][pixelIdx] *
           image[stampColIdx + xi + rPixX * (stampRowIdx + yi)];
    }
  }
  stamp->scprod[ncomp1 + 1] = q;

  return;
}

/**
 * @brief Perform an independent per-stamp kernel fit, sigma-clip outlier stamps
 *        by kernel sum, optionally do a global fit, and return a figure of
 * merit.
 *
 * @details This function is used to evaluate candidate convolution directions
 * ("which image is sharper?").  For each stamp it extracts the per-stamp
 * normal-equation system (stamp->mat, stamp->scprod), solves it independently
 * via LU decomposition, and records the kernel sum (the integral of the kernel,
 * equal to the flux-scaling ratio).  sigma_clip() is then applied to the
 * distribution of kernel sums; stamps deviating by more than kerSigReject sigma
 * are flagged.
 *
 * When forceConvolve == "b" (both directions tested) a global fit is performed
 * using only the unflagged stamps, and three figures of merit are computed
 * for each stamp:
 *   - merit1: mean chi-squared per pixel (variance metric)
 *   - merit2: standard deviation of the difference-image pixel distribution
 *   - merit3: histogram-FWHM-based noise estimate
 * All three are normalised by the kernel sum (so units are fractional noise
 * relative to the source flux).  The function returns the user-selected metric
 * (figMerit flag).
 *
 * @param stamps  Array of nS stamps with pre-built convolution vectors and
 *                per-stamp normal-equation data.
 * @param nS      Number of stamps.
 * @param imRef   Full-frame reference image.
 * @param imNoise Per-pixel noise image for chi-squared evaluation.
 * @return Sigma-normalised figure of merit (smaller is better); 0 if the
 *         global fit is not requested.
 */
double check_stamps(stamp_struct* stamps, int nS, float* imRef,
                    float* imNoise) {
  int nComps, stampIdx, matrixRowIdx, matrixColIdx, meritCount1, meritCount2,
      meritCount3;
  double sum = 0, kmean, kstdev;
  double merit1, merit2, merit3, sig1, sig2, sig3;
  float *meritValues1, *meritValues2, *meritValues3, *ks;
  int substampCenterX, substampCenterY, nks;

  double** matrix;
  int mat_size;
  int ncomp1, ncomp2, ncomp, nbg_vec;

  int ntestStamps;
  double* testKerSol = NULL;
  stamp_struct* testStamps = NULL;

  /* kernel sum */
  ks = (float*)calloc(nS, sizeof(float));
  nks = 0;

  ncomp1 = nCompKer - 1;
  ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  ncomp = ncomp1 * ncomp2;
  nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
  mat_size = ncomp1 * ncomp2 + nbg_vec + 1;

  LOG_DEBUG(" Mat_size0: %i ncomp2: %i ncomp1: %i nbg_vec: %i", mat_size,
            ncomp2, ncomp1, nbg_vec);

  /* for inital fit */
  nComps = nCompKer + 1;

  for (stampIdx = 0; stampIdx < nS; stampIdx++) {
    substampCenterX = stamps[stampIdx].xss[stamps[stampIdx].sscnt];
    substampCenterY = stamps[stampIdx].yss[stamps[stampIdx].sscnt];

    /* extract check_mat to solve one particular stamp */
    for (matrixRowIdx = 1; matrixRowIdx <= nComps; matrixRowIdx++) {
      check_vec[matrixRowIdx] = stamps[stampIdx].scprod[matrixRowIdx];

      for (matrixColIdx = 1; matrixColIdx <= matrixRowIdx; matrixColIdx++) {
        check_mat[matrixRowIdx][matrixColIdx] =
            stamps[stampIdx].mat[matrixRowIdx][matrixColIdx];
        check_mat[matrixColIdx][matrixRowIdx] =
            check_mat[matrixRowIdx][matrixColIdx];
      }
    }

    /* fit stamp, the constant kernel coefficients end up in check_vec */
    solve_spd(check_mat, nComps, check_vec);

    /* find kernel sum */
    sum = check_vec[1];
    check_stack[stampIdx] = sum;
    stamps[stampIdx].norm = sum;
    ks[nks++] = sum;

    LOG_DEBUG_INDENT(4, "# %d    xss: %4i yss: %4i  ksum: %f", stampIdx,
                     stamps[stampIdx].xss[stamps[stampIdx].sscnt],
                     stamps[stampIdx].yss[stamps[stampIdx].sscnt], sum);
  }

  sigma_clip(ks, nks, &kmean, &kstdev, 10);

  LOG_DEBUG_INDENT(4,
                   "%.1f sigma clipped mean ksum : %.3f, stdev : %.3f, n : %i",
                   kerSigReject, kmean, kstdev, nks);

  /* so we need some way to reject bad stamps here in the first test,
     we decided to use kernel sum.  is there a better way?  part of
     the trick is that if some things are variable, you get different
     kernel sums, but the subtraction itself should come out ok. */

  /* stamps.diff : delta ksum in sigma */

  /* here we want to reject high sigma points on the HIGH and LOW
     side, since we want things with the same normalization */
  for (stampIdx = 0; stampIdx < nS; stampIdx++) {
    stamps[stampIdx].diff = fabs((stamps[stampIdx].norm - kmean) / kstdev);
  }

  /*****************************************************
   * Global fit for kernel solution
   *****************************************************/

  /* do only if necessary */
  if ((strncmp(forceConvolve, "b", 1) == 0)) {
    /* allocate fitting matrix */
    matrix = (double**)calloc((mat_size + 1), sizeof(double*));
    for (stampIdx = 0; stampIdx <= mat_size; stampIdx++)
      matrix[stampIdx] = (double*)calloc((mat_size + 1), sizeof(double));

    /* allocate weight matrix */
    wxy = (double**)calloc(nS, sizeof(double*));
    for (stampIdx = 0; stampIdx < nS; stampIdx++)
      wxy[stampIdx] = (double*)calloc(ncomp2, sizeof(double));

    /* first find out how many good stamps to allocate */
    ntestStamps = 0;
    for (stampIdx = 0; stampIdx < nS; stampIdx++)
      if (stamps[stampIdx].diff < kerSigReject) {
        ntestStamps++;
      } else {
        LOG_DEBUG_INDENT(
            4, "# %d    skipping xss: %4i yss: %4i ksum: %f sigma: %f",
            stampIdx, stamps[stampIdx].xss[stamps[stampIdx].sscnt],
            stamps[stampIdx].yss[stamps[stampIdx].sscnt], stamps[stampIdx].norm,
            stamps[stampIdx].diff);
      }

    /* then allocate test stamp structure */
    if (!(testStamps =
              (stamp_struct*)calloc(ntestStamps, sizeof(stamp_struct)))) {
      LOG_ERROR("Failed to allocate test stamp list");
      exit(1);
    }
    testKerSol = (double*)calloc((nCompTotal + 1), sizeof(double));

    /* and point test stamp structure to good stamps */
    ntestStamps = 0;
    for (stampIdx = 0; stampIdx < nS; stampIdx++)
      if (stamps[stampIdx].diff < kerSigReject)
        testStamps[ntestStamps++] = stamps[stampIdx];

    /* finally do fit */
    LOG_DEBUG("Expanding Test Matrix For Fit");
    build_matrix(testStamps, ntestStamps, matrix);
    build_scprod(testStamps, ntestStamps, imRef, testKerSol);
    solve_spd(matrix, mat_size, testKerSol);

    /* get the kernel sum to normalize figures of merit! */
    kmean = make_kernel_dispatch(0, 0, testKerSol);

    /* determine figure of merit from good stamps */

    /* average of sum (diff**2 / value), ~variance */
    meritValues1 = (float*)calloc(ntestStamps, sizeof(float));

    /* standard deviation of pixel distribution */
    meritValues2 = (float*)calloc(ntestStamps, sizeof(float));

    /* noise sd based on histogram distribution width */
    meritValues3 = (float*)calloc(ntestStamps, sizeof(float));

    meritCount1 = 0;
    meritCount2 = 0;
    meritCount3 = 0;
    for (stampIdx = 0; stampIdx < ntestStamps; stampIdx++) {
      getStampSig(&testStamps[stampIdx], testKerSol, imNoise, &sig1, &sig2,
                  &sig3);

      if ((sig1 != -1) && (sig1 <= MAXVAL)) {
        meritValues1[meritCount1++] = sig1;
      }
      if ((sig2 != -1) && (sig2 <= MAXVAL)) {
        meritValues2[meritCount2++] = sig2;
      }
      if ((sig3 != -1) && (sig3 <= MAXVAL)) {
        meritValues3[meritCount3++] = sig3;
      }
    }
    sigma_clip(meritValues1, meritCount1, &merit1, &sig1, 10);
    sigma_clip(meritValues2, meritCount2, &merit2, &sig2, 10);
    sigma_clip(meritValues3, meritCount3, &merit3, &sig3, 10);

    /* normalize by kernel sum */
    merit1 /= kmean;
    merit2 /= kmean;
    merit3 /= kmean;

    /* clean up this mess */
    if (testKerSol) free(testKerSol);
    if (testStamps) free(testStamps);
    for (stampIdx = 0; stampIdx <= mat_size; stampIdx++) free(matrix[stampIdx]);
    for (stampIdx = 0; stampIdx < nS; stampIdx++) free(wxy[stampIdx]);
    free(matrix);
    free(wxy);

    free(meritValues1);
    free(meritValues2);
    free(meritValues3);
    free(ks);

    /* average value of figures of merit across stamps */
    LOG_PROGRESS("<var_merit> = %.3f, <sd_merit> = %.3f, <hist_merit> = %.3f",
                 merit1, merit2, merit3);

    /* return what is asked for if possible, if not use backup */
    if (strncmp(figMerit, "v", 1) == 0) {
      if (meritCount1 > 0) {
        return merit1;
      } else if (meritCount2 > 0) {
        return merit2;
      } else if (meritCount3 > 0) {
        return merit3;
      } else {
        return 666;
      }
    } else if (strncmp(figMerit, "s", 1) == 0) {
      if (meritCount2 > 0) {
        return merit2;
      } else if (meritCount1 > 0) {
        return merit1;
      } else if (meritCount3 > 0) {
        return merit3;
      } else {
        return 666;
      }
    } else if (strncmp(figMerit, "h", 1) == 0) {
      if (meritCount3 > 0) {
        return merit3;
      } else if (meritCount1 > 0) {
        return merit1;
      } else if (meritCount2 > 0) {
        return merit2;
      } else {
        return 666;
      }
    }
  } else
    return 0;

  return 0;
}

/**
 * @brief Assemble the global normal-equations matrix M for the
 * spatially-varying kernel fit by accumulating spatially-weighted per-stamp
 * cross-product matrices.
 *
 * @details Implements Alard & Lupton (1998), Eq. 6.  For each valid stamp
 * (those whose substamp counter has not overflowed), the spatial polynomial
 * weight vector wxy[istamp][k] = fx^i · fy^j is evaluated at the substamp
 * centre (normalised to [-1, 1] across the image).  The global matrix block
 * for the spatially-varying kernel coefficients is then:
 *   M[I][J] += wxy[i2] · wxy[j2] · stamp->mat[i1][j1]
 * where the composite indices I = (i1, i2) and J = (j1, j2) separate the
 * kernel basis index (i1) from the spatial polynomial index (i2).  Background
 * polynomial coupling blocks are accumulated similarly.  The full matrix is
 * symmetrised by copying the upper triangle to the lower triangle before
 * returning.
 *
 * @param stamps  Array of nS stamps with pre-built per-stamp normal-equation
 *                matrices (stamp->mat) and scalar products (stamp->scprod).
 * @param nS      Number of stamps.
 * @param matrix  Pre-allocated (mat_size+1)×(mat_size+1) output matrix;
 *                zeroed on entry and filled on return.
 */
void build_matrix(stamp_struct* stamps, int nS, double** matrix) {
  int mat_size, matrixRowIdx, matrixColIdx, pixStamp, stampIdx, polyTermIdx,
      gaussianCompIdx1, polyTermWithinGaussianIdx1, gaussianCompIdx2,
      polyTermIdx2, polyTermWithinGaussianIdx2, bgTermIdx, bgTermColIdx,
      bgVectorIdx, matrixColCalcIdx;
  int polyTermIdx1;
  int ncomp1, ncomp2, ncomp, nbg_vec;
  double **matrix0, p0, q;
  double** vec;
  float rPixX2, rPixY2;

  int polyDegX, polyDegY, stampCenterX, stampCenterY;
  double polyBasisX, polyBasisY, normalizedX, normalizedY;

  ncomp1 = nCompKer - 1;
  ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  ncomp = ncomp1 * ncomp2;
  nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;

  pixStamp = fwKSStamp * fwKSStamp;
  rPixX2 = 0.5 * rPixX;
  rPixY2 = 0.5 * rPixY;

  mat_size = ncomp1 * ncomp2 + nbg_vec + 1;
  LOG_DEBUG(" Mat_size: %i ncomp2: %i ncomp1: %i nbg_vec: %i", mat_size, ncomp2,
            ncomp1, nbg_vec);

  for (matrixRowIdx = 0; matrixRowIdx <= mat_size; matrixRowIdx++)
    for (matrixColIdx = 0; matrixColIdx <= mat_size; matrixColIdx++)
      matrix[matrixRowIdx][matrixColIdx] = 0.0;

  for (stampIdx = 0; stampIdx < nS; stampIdx++)
    for (polyTermIdx = 0; polyTermIdx < ncomp2; polyTermIdx++)
      wxy[stampIdx][polyTermIdx] = 0.0;

  for (stampIdx = 0; stampIdx < nS; stampIdx++) {
    /* skip over any bad stamps along the way */
    while (stamps[stampIdx].sscnt >= stamps[stampIdx].nss) {
      ++stampIdx;
      if (stampIdx >= nS) break;
    }
    if (stampIdx >= nS) break;

    vec = stamps[stampIdx].vectors;
    stampCenterX = stamps[stampIdx].xss[stamps[stampIdx].sscnt];
    stampCenterY = stamps[stampIdx].yss[stamps[stampIdx].sscnt];
    /* RANGE FROM -1 to 1 */
    normalizedX = (stampCenterX - rPixX2) / rPixX2;
    normalizedY = (stampCenterY - rPixY2) / rPixY2;

    /* build weight function *HERE* */
    polyTermIdx = 0;
    polyBasisX = 1.0;
    for (polyDegX = 0; polyDegX <= kerOrder; polyDegX++) {
      polyBasisY = 1.0;
      for (polyDegY = 0; polyDegY <= kerOrder - polyDegX; polyDegY++) {
        wxy[stampIdx][polyTermIdx++] = polyBasisX * polyBasisY;
        polyBasisY *= normalizedY;
      }
      polyBasisX *= normalizedX;
    }

    matrix0 = stamps[stampIdx].mat;
    for (polyTermIdx1 = 0; polyTermIdx1 < ncomp; polyTermIdx1++) {
      /* Decompose linear index polyTermIdx1 into (gaussianCompIdx1,
         polyTermWithinGaussianIdx1): gaussianCompIdx1 = which Gaussian
         component (0 to ncomp1-1) polyTermWithinGaussianIdx1 = which polynomial
         term within that component (0 to ncomp2-1) */
      gaussianCompIdx1 = polyTermIdx1 / ncomp2;
      polyTermWithinGaussianIdx1 = polyTermIdx1 - gaussianCompIdx1 * ncomp2;

      for (polyTermIdx2 = 0; polyTermIdx2 <= polyTermIdx1; polyTermIdx2++) {
        /* Same decomposition for column index polyTermIdx2 */
        gaussianCompIdx2 = polyTermIdx2 / ncomp2;
        polyTermWithinGaussianIdx2 = polyTermIdx2 - gaussianCompIdx2 * ncomp2;

        /* spatially weighted W_m1 and W_m2 integrals */
        matrix[(gaussianCompIdx1 * ncomp2 + polyTermWithinGaussianIdx1) + 2]
              [(gaussianCompIdx2 * ncomp2 + polyTermWithinGaussianIdx2) + 2] +=
            wxy[stampIdx][polyTermWithinGaussianIdx1] *
            wxy[stampIdx][polyTermWithinGaussianIdx2] *
            matrix0[gaussianCompIdx1 + 2][gaussianCompIdx2 + 2];
      }
    }

    matrix[1][1] += matrix0[1][1];
    for (polyTermIdx1 = 0; polyTermIdx1 < ncomp; polyTermIdx1++) {
      /* Index decomposition for background terms */
      gaussianCompIdx1 = polyTermIdx1 / ncomp2;
      polyTermWithinGaussianIdx1 = polyTermIdx1 - gaussianCompIdx1 * ncomp2;
      matrix[(gaussianCompIdx1 * ncomp2 + polyTermWithinGaussianIdx1) + 2][1] +=
          wxy[stampIdx][polyTermWithinGaussianIdx1] *
          matrix0[gaussianCompIdx1 + 2][1];
    }

    for (bgTermIdx = 0; bgTermIdx < nbg_vec; bgTermIdx++) {
      polyTermIdx1 = ncomp + bgTermIdx + 1;
      bgVectorIdx = ncomp1 + bgTermIdx + 1;
      for (gaussianCompIdx1 = 1; gaussianCompIdx1 < ncomp1 + 1;
           gaussianCompIdx1++) {
        p0 = 0.0;

        /* integrate convolved images over all order backgrounds */
        for (polyTermIdx = 0; polyTermIdx < pixStamp; polyTermIdx++)
          p0 += vec[gaussianCompIdx1][polyTermIdx] *
                vec[bgVectorIdx][polyTermIdx];

        /* spatially weighted image * background terms */
        for (polyTermIdx2 = 0; polyTermIdx2 < ncomp2; polyTermIdx2++) {
          matrixColCalcIdx = (gaussianCompIdx1 - 1) * ncomp2 + polyTermIdx2 + 1;
          matrix[polyTermIdx1 + 1][matrixColCalcIdx + 1] +=
              p0 * wxy[stampIdx][polyTermIdx2];
        }
      }

      p0 = 0.0;
      for (polyTermIdx = 0; polyTermIdx < pixStamp; polyTermIdx++)
        p0 += vec[0][polyTermIdx] * vec[bgVectorIdx][polyTermIdx];
      matrix[polyTermIdx1 + 1][1] += p0;

      /* background * background */
      for (bgTermColIdx = 0; bgTermColIdx <= bgTermIdx; bgTermColIdx++) {
        for (polyTermIdx = 0, q = 0.0; polyTermIdx < pixStamp; polyTermIdx++)
          q += vec[bgVectorIdx][polyTermIdx] *
               vec[ncomp1 + bgTermColIdx + 1][polyTermIdx];
        matrix[polyTermIdx1 + 1][ncomp + bgTermColIdx + 2] += q;
      }
    }
  }

  /* fill lower half of matrix */
  for (matrixRowIdx = 0; matrixRowIdx < mat_size; matrixRowIdx++) {
    for (matrixColIdx = 0; matrixColIdx <= matrixRowIdx; matrixColIdx++) {
      matrix[matrixColIdx + 1][matrixRowIdx + 1] =
          matrix[matrixRowIdx + 1][matrixColIdx + 1];
      /* fprintf(stderr, "matrix[%i][%i]: %lf\n", i,j,matrix[i+1][j+1]); */
    }
  }

  return;
}

/**
 * @brief Assemble the global right-hand-side vector b for the normal equations
 *        by accumulating spatially-weighted per-stamp scalar products.
 *
 * @details For each valid stamp the pre-computed scalar products
 * stamp->scprod[i] (built by build_scprod0()) are scaled by the appropriate
 * spatial polynomial weight wxy[istamp][i2] and accumulated into @p kernelSol:
 *   kernelSol[I] += wxy[istamp][i2] · stamp->scprod[i1]
 * Background terms are computed directly from the background basis images and
 * the reference pixel values.  On entry @p kernelSol is zeroed.  After
 * build_matrix() fills M and this function fills b, LAPACK Cholesky solves
 * M·a=b to yield the kernel coefficient vector in @p kernelSol.
 *
 * @param stamps     Array of nS stamps.
 * @param nS         Number of stamps.
 * @param image      Full-frame reference image (used for the background terms).
 * @param kernelSol  Output: right-hand-side vector of length nCompTotal+1;
 *                   zeroed on entry, filled on return.
 */
void build_scprod(stamp_struct* stamps, int nS, float* image,
                  double* kernelSol) {
  int stampIdx, stampColIdx, stampRowIdx, stampCenterX, stampCenterY,
      gaussianCompIdx, polyTermIdx, pixelIdx, bgVecIdx, solIdx,
      recomputedColIdx;
  int ncomp1, ncomp2, ncomp, nbg_vec;
  double scalarProduct, innerProduct;
  double** vectorArray;

  ncomp1 = nCompKer - 1;
  ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  ncomp = ncomp1 * ncomp2;
  nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;

  for (solIdx = 0; solIdx <= ncomp + nbg_vec + 1; solIdx++)
    kernelSol[solIdx] = 0.0;

  for (stampIdx = 0; stampIdx < nS; stampIdx++) {
    /* skip over any bad stamps along the way */
    while (stamps[stampIdx].sscnt >= stamps[stampIdx].nss) {
      ++stampIdx;
      if (stampIdx >= nS) break;
    }
    if (stampIdx >= nS) break;

    vectorArray = stamps[stampIdx].vectors;
    stampCenterX = stamps[stampIdx].xss[stamps[stampIdx].sscnt];
    stampCenterY = stamps[stampIdx].yss[stamps[stampIdx].sscnt];

    scalarProduct = stamps[stampIdx].scprod[1];
    kernelSol[1] += scalarProduct;

    /* spatially weighted convolved image * ref image */
    for (gaussianCompIdx = 1; gaussianCompIdx < ncomp1 + 1; gaussianCompIdx++) {
      scalarProduct = stamps[stampIdx].scprod[gaussianCompIdx + 1];
      for (polyTermIdx = 0; polyTermIdx < ncomp2; polyTermIdx++) {
        recomputedColIdx = (gaussianCompIdx - 1) * ncomp2 + polyTermIdx + 1;
        /* no need for weighting here */
        kernelSol[recomputedColIdx + 1] +=
            scalarProduct * wxy[stampIdx][polyTermIdx];
      }
    }

    /* spatially weighted bg model convolved with ref image */
    for (bgVecIdx = 0; bgVecIdx < nbg_vec; bgVecIdx++) {
      innerProduct = 0.0;
      for (stampColIdx = -hwKSStamp; stampColIdx <= hwKSStamp; stampColIdx++) {
        for (stampRowIdx = -hwKSStamp; stampRowIdx <= hwKSStamp;
             stampRowIdx++) {
          pixelIdx =
              stampColIdx + hwKSStamp + fwKSStamp * (stampRowIdx + hwKSStamp);
          innerProduct += vectorArray[ncomp1 + bgVecIdx + 1][pixelIdx] *
                          image[stampColIdx + stampCenterX +
                                rPixX * (stampRowIdx + stampCenterY)];
        }
      }
      kernelSol[ncomp + bgVecIdx + 2] += innerProduct;
    }
  }
  return;
}

/**
 * @brief Compute the mean chi-squared per pixel for a substamp in the final
 *        difference image.
 *
 * @details Sums (imDiff[k])² / (imNoise[k])² over all unmasked pixels in the
 * fwKSStamp×fwKSStamp substamp centred on the stamp's active substamp location
 * and divides by the pixel count to produce a mean chi-squared.  Pixels
 * flagged as FLAG_INPUT_ISBAD are skipped.  Returns -1 via @p sig if no valid
 * pixels are found.
 *
 * @param stamp    Stamp whose active substamp (stamp->sscnt) defines the
 *                 evaluation region.
 * @param imDiff   Full-frame difference image.
 * @param imNoise  Full-frame per-pixel noise image (sigma values).
 * @param sig      Output: mean chi-squared per pixel over the substamp, or -1
 *                 if no valid pixels exist.
 */
void getFinalStampSig(stamp_struct* stamp, float* imDiff, float* imNoise,
                      double* sig) {
  int colIdx, rowIdx, pixelIdx, validPixelCount = 0;
  int xPixel, xRegion = stamp->xss[stamp->sscnt];
  int yPixel, yRegion = stamp->yss[stamp->sscnt];
  float diffData, noiseInv;

  *sig = 0;

  for (rowIdx = 0; rowIdx < fwKSStamp; rowIdx++) {
    yPixel = yRegion - hwKSStamp + rowIdx;

    for (colIdx = 0; colIdx < fwKSStamp; colIdx++) {
      xPixel = xRegion - hwKSStamp + colIdx;

      pixelIdx = xPixel + rPixX * yPixel;
      diffData = imDiff[pixelIdx];
      noiseInv = 1. / imNoise[pixelIdx];

      /* this shouldn't be the case, but just in case... */
      if (mRData[pixelIdx] & FLAG_INPUT_ISBAD) continue;

      validPixelCount++;
      *sig += diffData * diffData * noiseInv * noiseInv;
    }
  }
  if (validPixelCount > 0)
    *sig /= validPixelCount;
  else
    *sig = -1;

  return;
}

/**
 * @brief Evaluate three quality-of-fit metrics for a single substamp given the
 *        current kernel solution.
 *
 * @details Builds the model convolved image for the substamp via make_model()
 * and subtracts the reference (stamp->krefArea), then computes:
 *  - @p sig1: mean chi-squared per pixel (diff²/noise²), normalised by pixel
 *    count; used as the "variance" figure of merit (figMerit=="v").
 *  - @p sig2: standard deviation of the residual pixel distribution from
 *    getStampStats3(); used as the "sd" figure of merit (figMerit=="s").
 *  - @p sig3: histogram-FWHM noise estimate from getStampStats3(); used as the
 *    "histogram" figure of merit (figMerit=="h").
 * Any metric that cannot be computed (e.g. masked pixels, NaN values) is
 * returned as -1.
 *
 * @param stamp      Stamp to evaluate; stamp->krefArea must be filled.
 * @param kernelSol  Current global kernel solution vector.
 * @param imNoise    Full-frame per-pixel noise image.
 * @param sig1       Output: mean chi-squared per pixel metric.
 * @param sig2       Output: pixel-distribution standard deviation metric.
 * @param sig3       Output: histogram FWHM metric.
 */
void getStampSig(stamp_struct* stamp, double* kernelSol, float* imNoise,
                 double* sig1, double* sig2, double* sig3) {
  int colIdx, rowIdx, stampPixelIdx, substampIdx, validPixelCount, xRegion,
      yRegion, xPixel, yPixel;
  double cSum, cMean, cMedian, cMode, cLfwhm;
  double *referenceImage, templateData, imageData, noiseData, diffValue,
      bgValue;

  /* info */
  substampIdx = stamp->sscnt;
  xRegion = stamp->xss[substampIdx];
  yRegion = stamp->yss[substampIdx];

  /* the comparison image */
  referenceImage = stamp->krefArea;
  /* background from fit */
  bgValue = get_background(xRegion, yRegion, kernelSol);
  /* temp contains the convolved image from fit, fwKSStamp x fwKSStamp */
  make_model(stamp, kernelSol, temp);

  /* get sigma of stamp diff */
  validPixelCount = 0;
  *sig1 = 0;
  *sig2 = 0;
  *sig3 = 0;

  for (rowIdx = 0; rowIdx < fwKSStamp; rowIdx++) {
    yPixel = yRegion - hwKSStamp + rowIdx;

    for (colIdx = 0; colIdx < fwKSStamp; colIdx++) {
      xPixel = xRegion - hwKSStamp + colIdx;

      stampPixelIdx = colIdx + rowIdx * fwKSStamp;

      templateData = temp[stampPixelIdx];
      imageData = referenceImage[stampPixelIdx];
      noiseData = imNoise[xPixel + rPixX * yPixel];

      diffValue = templateData - imageData + bgValue;

      if ((mRData[xPixel + rPixX * yPixel] & FLAG_INPUT_ISBAD) ||
          (fabs(imageData) <= ZEROVAL)) {
        continue;
      } else {
        temp[stampPixelIdx] = diffValue;
      }

      /* check for NaN */
      if ((templateData * 0.0 != 0.0) || (imageData * 0.0 != 0.0)) {
        mRData[xPixel + rPixX * yPixel] |= (FLAG_INPUT_ISBAD | FLAG_ISNAN);
        continue;
      }

      validPixelCount++;
      *sig1 += diffValue * diffValue / noiseData;
      /*fprintf(stderr, "OK %d %d : %f %f %f\n", xPixel, yPixel, templateData,
       * imageData, noiseData);*/
    }
  }
  if (validPixelCount > 0) {
    *sig1 /= validPixelCount;
    if (*sig1 >= MAXVAL) *sig1 = -1;
  } else
    *sig1 = -1;

  /* don't do think unless you need to! */
  if (strncmp(figMerit, "v", 1) != 0) {
    if (getStampStats3(temp, xRegion - hwKSStamp, yRegion - hwKSStamp,
                       fwKSStamp, fwKSStamp, &cSum, &cMean, &cMedian, &cMode,
                       sig2, sig3, &cLfwhm, 0x0, 0xffff, 5)) {
      *sig2 = -1;
      *sig3 = -1;
    } else if (*sig2 < 0 || *sig2 >= MAXVAL)
      *sig2 = -1;
    else if (*sig3 < 0 || *sig3 >= MAXVAL)
      *sig3 = -1;
  }

  return;
}

/**
 * @brief Sigma-clip stamps after a global kernel fit and flag those with poor
 *        residuals for substamp replacement.
 *
 * @details For each stamp that is currently represented by a valid substamp,
 * getStampSig() is called to measure the residual quality metric selected by
 * the figMerit flag.  Stamps with an invalid metric (sig == -1) are immediately
 * advanced to the next substamp via fillStamp().  Valid metrics are collected
 * and sigma-clipped to derive a robust mean and scatter; any stamp whose metric
 * exceeds mean + kerSigReject·stdev is likewise advanced to its next substamp.
 * If any stamp was advanced the function returns 1 (triggering a matrix rebuild
 * and re-solve in fitKernel()); otherwise it returns 0 (convergence).
 *
 * The sigma-clipped mean and scatter are stored in @p meansigSubstamps and
 * @p scatterSubstamps for inclusion in the output FITS header.
 *
 * @param stamps              Array of nS stamps.
 * @param kernelSol           Current kernel solution vector.
 * @param imConv              Full-frame image to be convolved (used by
 *                            fillStamp() when advancing to a new substamp).
 * @param imRef               Full-frame reference image.
 * @param imNoise             Full-frame noise image for chi-squared evaluation.
 * @param meansigSubstamps    Output: sigma-clipped mean residual at the end of
 *                            this iteration.
 * @param scatterSubstamps    Output: sigma-clipped scatter of residuals.
 * @param NskippedSubstamps   Output: number of stamps skipped (out of
 *                            substamps).
 * @return 1 if any stamp was advanced (another iteration needed), 0 if
 *         converged.
 */
char check_again(stamp_struct* stamps, double* kernelSol, float* imConv,
                 float* imRef, float* imNoise, double* meansigSubstamps,
                 double* scatterSubstamps, int* NskippedSubstamps) {
  int stampIdx, validStampCount, goodStampCount;
  double qualityMetric, mean, stdev;
  char needsReiteration;
  double sig1, sig2, sig3;
  float* qualityMetrics;

  qualityMetrics = (float*)calloc(nS, sizeof(float));
  validStampCount = 0;

  qualityMetric = 0;
  needsReiteration = 0;
  mean = stdev = 0.0;
  *NskippedSubstamps = 0;

  for (stampIdx = 0; stampIdx < nS; stampIdx++) {
    /* if was fit with a good legit substamp */
    if (stamps[stampIdx].sscnt < stamps[stampIdx].nss) {
      getStampSig(&stamps[stampIdx], kernelSol, imNoise, &sig1, &sig2, &sig3);

      if ((strncmp(figMerit, "v", 1) == 0 && (sig1 == -1)) ||
          (strncmp(figMerit, "s", 1) == 0 && (sig2 == -1)) ||
          (strncmp(figMerit, "h", 1) == 0 && (sig3 == -1))) {
        /* something went wrong with this one... */
        LOG_DEBUG_INDENT(4,
                         "# %d    xss: %4i yss: %4i sig: %6.3f sscnt: %2i nss: "
                         "%2i ITERATE substamp (BAD)",
                         stampIdx, stamps[stampIdx].xss[stamps[stampIdx].sscnt],
                         stamps[stampIdx].yss[stamps[stampIdx].sscnt],
                         qualityMetric, stamps[stampIdx].sscnt,
                         stamps[stampIdx].nss);

        stamps[stampIdx].sscnt++;
        fillStamp(&stamps[stampIdx], imConv, imRef);

        needsReiteration = 1;

      } else {
        if (strncmp(figMerit, "v", 1) == 0)
          qualityMetric = sig1;
        else if (strncmp(figMerit, "s", 1) == 0)
          qualityMetric = sig2;
        else if (strncmp(figMerit, "h", 1) == 0)
          qualityMetric = sig3;

        LOG_DEBUG_INDENT(
            4, "# %d    xss: %4i yss: %4i sig: %6.3f sscnt: %2i nss: %2i OK",
            stampIdx, stamps[stampIdx].xss[stamps[stampIdx].sscnt],
            stamps[stampIdx].yss[stamps[stampIdx].sscnt], qualityMetric,
            stamps[stampIdx].sscnt, stamps[stampIdx].nss);

        stamps[stampIdx].chi2 = qualityMetric;
        qualityMetrics[validStampCount++] = qualityMetric;
      }
    } else {
      (*NskippedSubstamps)++;
      LOG_DEBUG_INDENT(4, "xs : %4i ys : %4i skipping...", stamps[stampIdx].x,
                       stamps[stampIdx].y);
    }
  }

  sigma_clip(qualityMetrics, validStampCount, &mean, &stdev, 10);
  LOG_PROGRESS("Mean sig: %6.3f stdev: %6.3f", mean, stdev);
  LOG_PROGRESS("Iterating through stamps with sig > %.3f",
               mean + kerSigReject * stdev);

  /* save the mean and scatter so that it can be saved in the fits header */
  (*meansigSubstamps) = mean;
  (*scatterSubstamps) = stdev;

  goodStampCount = 0;
  for (stampIdx = 0; stampIdx < nS; stampIdx++) {
    /* if currently represented by a good substamp */
    if (stamps[stampIdx].sscnt < stamps[stampIdx].nss) {
      /* no fabs() here, keep good stamps kerSigReject on the low side! */
      if ((stamps[stampIdx].chi2 - mean) > kerSigReject * stdev) {
        LOG_DEBUG_INDENT(4,
                         "# %d    xss: %4i yss: %4i sig: %6.3f sscnt: %2i nss: "
                         "%2i ITERATE substamp (poor sig)",
                         stampIdx, stamps[stampIdx].xss[stamps[stampIdx].sscnt],
                         stamps[stampIdx].yss[stamps[stampIdx].sscnt],
                         stamps[stampIdx].chi2, stamps[stampIdx].sscnt,
                         stamps[stampIdx].nss);

        stamps[stampIdx].sscnt++;
        goodStampCount += (!(fillStamp(&stamps[stampIdx], imConv, imRef)));

        needsReiteration = 1;
      } else
        goodStampCount += 1;
    }
  }

  LOG_PROGRESS("%d out of %d stamps remain", goodStampCount, nS);

  free(qualityMetrics);
  return needsReiteration;
}

/**
 * @brief Thread-safe kernel evaluator writing into caller-supplied buffers.
 *
 * @details Identical logic to make_kernel() but writes to lkernel and
 * lkernel_coeffs instead of the global arrays, allowing concurrent calls from
 * multiple OpenMP threads each with their own private buffers.
 *
 * @param xi             x pixel coordinate.
 * @param yi             y pixel coordinate.
 * @param kernelSol      Kernel solution vector (read-only).
 * @param lkernel        Caller-allocated buffer of size fwKernel*fwKernel.
 * @param lkernel_coeffs Caller-allocated buffer of size nCompKer.
 * @return Sum of all pixels in the assembled kernel image.
 */
static double make_kernel_local(int xi, int yi, double* kernelSol,
                                double* lkernel, double* lkernel_coeffs) {
  int gaussianCompIdx, solutionIdx, polyDegX, polyDegY, kernelPixelIdx;
  double polyBasisX, polyBasisY, kernelSum;
  double normalizedX, normalizedY;

  solutionIdx = 2;
  normalizedX = (xi - 0.5 * rPixX) / (0.5 * rPixX);
  normalizedY = (yi - 0.5 * rPixY) / (0.5 * rPixY);

  for (gaussianCompIdx = 1; gaussianCompIdx < nCompKer; gaussianCompIdx++) {
    lkernel_coeffs[gaussianCompIdx] = 0.0;
    polyBasisX = 1.0;
    for (polyDegX = 0; polyDegX <= kerOrder; polyDegX++) {
      polyBasisY = 1.0;
      for (polyDegY = 0; polyDegY <= kerOrder - polyDegX; polyDegY++) {
        lkernel_coeffs[gaussianCompIdx] +=
            kernelSol[solutionIdx++] * polyBasisX * polyBasisY;
        polyBasisY *= normalizedY;
      }
      polyBasisX *= normalizedX;
    }
  }
  lkernel_coeffs[0] = kernelSol[1];

  for (kernelPixelIdx = 0; kernelPixelIdx < fwKernel * fwKernel;
       kernelPixelIdx++)
    lkernel[kernelPixelIdx] = 0.0;

  kernelSum = 0.0;
  for (kernelPixelIdx = 0; kernelPixelIdx < fwKernel * fwKernel;
       kernelPixelIdx++) {
    for (gaussianCompIdx = 0; gaussianCompIdx < nCompKer; gaussianCompIdx++)
      lkernel[kernelPixelIdx] += lkernel_coeffs[gaussianCompIdx] *
                                 kernel_vec[gaussianCompIdx][kernelPixelIdx];
    kernelSum += lkernel[kernelPixelIdx];
  }
  return kernelSum;
}

/**
 * @brief TPS-based alternative to make_kernel_local().
 *
 * @details Evaluates TPS-interpolated kernel coefficients and assembles
 * the kernel into the local buffer (for FFT-accelerated convolution).
 *
 * @param xi              x pixel coordinate
 * @param yi              y pixel coordinate
 * @param kernelSol       Extended kernel solution (includes TPS parameters)
 * @param lkernel         Local kernel buffer (output, fwKernel×fwKernel)
 * @param lkernel_coeffs  Local coefficient buffer (output, nCompKer values)
 * @return Kernel sum (photometric scaling factor)
 */
static double make_kernel_local_tps(int xi, int yi, double* kernelSol,
                                    double* lkernel,
                                    double* lkernel_coeffs) {
  int gaussianCompIdx, kernelPixelIdx;
  double kernelSum;
  double *tps_weights, *tps_poly, *positions;

  /* Evaluate TPS-interpolated coefficients at (xi, yi) */
  for (gaussianCompIdx = 1; gaussianCompIdx < nCompKer; gaussianCompIdx++) {
    tps_weights = kernelSol + kernelSol_offset_tps_weights(gaussianCompIdx, nS);
    tps_poly = kernelSol + kernelSol_offset_tps_poly(gaussianCompIdx, nS);
    positions = kernelSol + kernelSol_offset_tps_positions(nS);

    lkernel_coeffs[gaussianCompIdx] =
        tps_evaluate((double)xi, (double)yi, positions, nS, tps_weights, tps_poly);
  }

  /* Constant component (no spatial variation) */
  lkernel_coeffs[0] = kernelSol[1];

  /* Assemble kernel from basis functions */
  for (kernelPixelIdx = 0; kernelPixelIdx < fwKernel * fwKernel;
       kernelPixelIdx++)
    lkernel[kernelPixelIdx] = 0.0;

  kernelSum = 0.0;
  for (kernelPixelIdx = 0; kernelPixelIdx < fwKernel * fwKernel;
       kernelPixelIdx++) {
    for (gaussianCompIdx = 0; gaussianCompIdx < nCompKer; gaussianCompIdx++)
      lkernel[kernelPixelIdx] += lkernel_coeffs[gaussianCompIdx] *
                                 kernel_vec[gaussianCompIdx][kernelPixelIdx];
    kernelSum += lkernel[kernelPixelIdx];
  }
  return kernelSum;
}

/**
 * @brief Dispatcher: evaluate kernel into local buffer using polynomial or TPS.
 *
 * @param xi              x pixel coordinate
 * @param yi              y pixel coordinate
 * @param kernelSol       Kernel solution vector
 * @param lkernel         Local kernel buffer (output)
 * @param lkernel_coeffs  Local coefficient buffer (output)
 * @return Kernel sum
 */
static double make_kernel_local_dispatch(int xi, int yi, double* kernelSol,
                                         double* lkernel,
                                         double* lkernel_coeffs) {
  if (useTPS) {
    return make_kernel_local_tps(xi, yi, kernelSol, lkernel, lkernel_coeffs);
  } else {
    return make_kernel_local(xi, yi, kernelSol, lkernel, lkernel_coeffs);
  }
}

/**
 * @brief Round n up to the smallest integer >= n whose only prime factors
 *        are 2, 3, or 5.  Such sizes have fast FFTW r2c/c2r plans.
 */
static int next_fftw_size(int n) {
  int m;
  for (;;) {
    m = n;
    while (m % 2 == 0) m /= 2;
    while (m % 3 == 0) m /= 3;
    while (m % 5 == 0) m /= 5;
    if (m == 1) return n;
    n++;
  }
}

/**
 * @brief FFT-accelerated implementation of spatial_convolve().
 *
 * @details Uses the linearity of convolution to restructure the spatially-
 * varying kernel as a small set of fixed-kernel convolutions followed by a
 * per-pixel polynomial-weighted sum (Alard & Lupton 1998, §2):
 *
 *   K(x,y) = Σᵢ cᵢ(x,y)·φᵢ
 *
 * Writing cᵢ(x,y) = Σⱼ aᵢⱼ·Pⱼ(x,y) and grouping by polynomial term j:
 *
 *   output(x,y) = kernelSol[1] · (I⊗φ₀)[x,y]
 *               + Σⱼ Pⱼ(x,y) · (I⊗effKernelⱼ)[x,y]
 *
 * where effKernelⱼ = Σᵢ kernelSol[2+(i-1)·ncomp2+j] · φᵢ.
 *
 * This requires only 1+ncomp2 FFT convolutions (7 for default kerOrder=2)
 * instead of nCompKer ≈ 50.  Memory cost: (1+ncomp2)·xSize·ySize·sizeof(float).
 * Time complexity: O(N log N) vs O(N·k²) for direct convolution (3-8× faster).
 *
 * Buffer Management (Critical for multi-region processing):
 *   - real_buf: Allocated as fftw_malloc(fft_nx * fft_ny * sizeof(double))
 *     Input/output buffer for FFTW r2c/c2r transforms. Padded to fft_nx×fft_ny
 *     to avoid circular wrap-around artifacts in periodic convolution.
 *
 *   - img_fft, ker_fft: Each allocated as fftw_malloc(nc_fft * sizeof(fftw_complex))
 *     where nc_fft = fft_ny * (fft_nx / 2 + 1) per FFTW 2D r2c output shape.
 *     CRITICAL: This formula MUST match plan shape (fft_ny, fft_nx).
 *     Bug history: Previous code used nc_fft = fft_nx * (fft_ny / 2 + 1),
 *     causing 20-element underallocation for non-square FFTs and heap corruption.
 *
 *   - effConv: Array of (1+ncomp2) float pointers, each pointing to xSize×ySize
 *     region (not padded). Stores the effective convolution results before
 *     polynomial weighting and assembly into cRdata.
 *
 * FFTW plan creation is serialised via #pragma omp critical(fftw_plan);
 * plan execution and the pixel-sum loop are parallelised with OpenMP.
 * Mask propagation uses direct loop (not FFT) to preserve binary mask semantics.
 * CFITSIO is not called here; all FITS I/O occurs in main.c.
 *
 * @param image      Full-frame input image (region pixel buffer, xSize×ySize).
 * @param variance   Pointer to variance image; updated in-place if non-NULL.
 * @param xSize      Region width in pixels (actual region, not padded).
 * @param ySize      Region height in pixels (actual region, not padded).
 * @param kernelSol  Kernel solution vector from fitKernel().
 * @param cRdata     Output convolved image (pre-allocated by caller, xSize×ySize).
 * @param cMask      Input pixel mask for mask propagation (xSize×ySize).
 */
static void spatial_convolve_fft(float* image, float** variance, int xSize,
                                 int ySize, double* kernelSol, float* cRdata,
                                 int* cMask) {
  /* --- variable declarations (all at top for safe goto use) --- */
  int ncomp2, nEffConv, fft_nx, fft_ny, nc_fft;
  int i, j, p, ik, jk, ii, jj, xp, yp, ni, ix, iy;
  int i1, i2, j2, i0, j0, j1, nsteps_x_loc, nsteps_y_loc;
  int ic, jc, nc, mbit, dovar, cidx;
  double ax, ay, out, xf, yf, kk, aks, uks, qv, norm, rPixX2, rPixY2;
  double *effKernel, *real_buf, *lkernel, *lkernel_coeffs;
  fftw_complex *img_fft, *ker_fft;
  fftw_plan plan_fwd, plan_inv;
  float **effConv, *vData;
  size_t image_size, effconv_element_size;

  effKernel = NULL;
  real_buf = NULL;
  lkernel = NULL;
  lkernel_coeffs = NULL;
  img_fft = NULL;
  ker_fft = NULL;
  plan_fwd = NULL;
  plan_inv = NULL;
  effConv = NULL;
  vData = NULL;

  /* Sanity checks on input dimensions */
  image_size = (size_t)xSize * ySize;
  if (xSize <= 0 || ySize <= 0 || image_size == 0) {
    LOG_ERROR("spatial_convolve_fft: invalid dimensions xSize=%d, ySize=%d",
              xSize, ySize);
    return;
  }
  if (!image || !cRdata || !cMask) {
    LOG_ERROR("spatial_convolve_fft: NULL input pointer");
    return;
  }

  ncomp2 = ((kerOrder + 1) * (kerOrder + 2)) / 2;
  nEffConv = 1 + ncomp2;
  dovar = (*variance != NULL);

  if (dovar) {
    vData = (float*)calloc((size_t)xSize * ySize, sizeof(float));
    if (!vData) {
      LOG_ERROR("Spatial_convolve_fft: out of memory (vData)");
      return;
    }
  }

  /* --- FFT dimensions: pad to avoid circular wrap-around --- */
  fft_nx = next_fftw_size(xSize + fwKernel - 1);
  fft_ny = next_fftw_size(ySize + fwKernel - 1);

  /* FFTW 2D Transform Documentation:
   *
   * FFTW 2D transforms are specified as fftw_plan_dft_r2c_2d(n0, n1, ...)
   * where (n0, n1) = (rows, cols) in C row-major order.
   *
   * In this code:
   *   - Shape parameter: (fft_ny, fft_nx) = (rows, cols)
   *   - Real input buffer: fft_ny rows × fft_nx cols (total: fft_ny*fft_nx doubles)
   *   - Complex output (r2c): fft_ny rows × (fft_nx/2+1) cols (total: fft_ny*(fft_nx/2+1) complexes)
   *   - Real output (c2r): fft_ny rows × fft_nx cols (total: fft_ny*fft_nx doubles, pre-allocated)
   *
   * Key Points:
   *   1. Memory layout is C row-major: element [i,j] is at index i*cols + j
   *   2. For r2c, output columns = (input_cols / 2) + 1 (Hermitian symmetry)
   *   3. For c2r, output size is same as input (padded) size
   *
   * Buffer Allocation:
   *   - real_buf:  fftw_malloc(fft_nx * fft_ny * sizeof(double))
   *   - img_fft:   fftw_malloc(nc_fft * sizeof(fftw_complex)) where
   *   - ker_fft:   nc_fft = fft_ny * (fft_nx / 2 + 1)
   *
   * CRITICAL: This formula MUST match the FFTW plan shape (fft_ny, fft_nx).
   * Previous bug: nc_fft = fft_nx * (fft_ny / 2 + 1) caused 20-element underallocation
   * for 280x300 FFT, leading to heap corruption in multi-region processing.
   */
  nc_fft = fft_ny * (fft_nx / 2 + 1);

  real_buf = (double*)fftw_malloc((size_t)fft_nx * fft_ny * sizeof(double));
  img_fft = (fftw_complex*)fftw_malloc((size_t)nc_fft * sizeof(fftw_complex));
  ker_fft = (fftw_complex*)fftw_malloc((size_t)nc_fft * sizeof(fftw_complex));
  effKernel = (double*)malloc((size_t)fwKernel * fwKernel * sizeof(double));

  if (!real_buf || !img_fft || !ker_fft || !effKernel) {
    LOG_ERROR("spatial_convolve_fft: out of memory (FFT buffers)");
    goto cleanup_fft;
  }

  effConv = (float**)calloc((size_t)nEffConv, sizeof(float*));
  if (!effConv) {
    LOG_ERROR("spatial_convolve_fft: out of memory (effConv)");
    goto cleanup_fft;
  }
  /* Initialize all pointers to NULL for safe cleanup if allocation fails */
  for (i = 0; i < nEffConv; i++) {
    effConv[i] = NULL;
  }
  /* Now allocate the actual buffers */
  for (i = 0; i < nEffConv; i++) {
    effConv[i] = (float*)malloc((size_t)xSize * ySize * sizeof(float));
    if (!effConv[i]) {
      LOG_ERROR("spatial_convolve_fft: out of memory (effConv[%d])", i);
      goto cleanup_fft;
    }
  }

  /* --- Create FFTW plans (not thread-safe: serialise with critical) ---
   * Use FFTW_MEASURE (optimized, slow planning) in single-region mode, or
   * FFTW_ESTIMATE (quick planning) in multi-region mode.
   */
#ifdef _OPENMP
#pragma omp critical(fftw_plan)
#endif
  {
    int plan_flags = (nRegX == 1 && nRegY == 1) ? FFTW_MEASURE : FFTW_ESTIMATE;
    plan_fwd = fftw_plan_dft_r2c_2d(fft_ny, fft_nx, real_buf, ker_fft, plan_flags);
    plan_inv = fftw_plan_dft_c2r_2d(fft_ny, fft_nx, ker_fft, real_buf, plan_flags);
  }

  if (!plan_fwd || !plan_inv) {
    LOG_ERROR("spatial_convolve_fft: FFTW plan creation failed");
    goto cleanup_fft;
  }

  /* --- FFT the image once; save spectrum in img_fft ---
   * Copy image data into padded real_buf, leaving the padding region as zeros.
   * Buffer indexing: real_buf and image are stored in C row-major order:
   *   element [y][x] is at index x + stride*y
   * Stride: real_buf has stride fft_nx (padded), image has stride xSize (unpadded).
   * We copy only the valid region [0:ySize)[0:xSize), padding with zeros.
   */
  memset(real_buf, 0, (size_t)fft_nx * fft_ny * sizeof(double));
  for (yp = 0; yp < ySize; yp++)
    for (xp = 0; xp < xSize; xp++)
      real_buf[xp + fft_nx * yp] = image[xp + xSize * yp];

  fftw_execute_dft_r2c(plan_fwd, real_buf, ker_fft);
  memcpy(img_fft, ker_fft, (size_t)nc_fft * sizeof(fftw_complex));

  /* --- Compute nEffConv effective basis convolutions ---
   * j == -1 : constant term  -> effConv[0]
   * j == 0..ncomp2-1 : poly terms -> effConv[1..ncomp2] */
  for (j = -1; j < ncomp2; j++) {
    cidx = j + 1;

    memset(effKernel, 0, (size_t)fwKernel * fwKernel * sizeof(double));
    if (j < 0) {
      /* constant term: kernelSol[1] * kernel_vec[0] */
      for (p = 0; p < fwKernel * fwKernel; p++)
        effKernel[p] = kernelSol[1] * kernel_vec[0][p];
    } else {
      /* poly term j: Σᵢ kernelSol[2+(i-1)*ncomp2+j] * kernel_vec[i] */
      for (i1 = 1; i1 < nCompKer; i1++) {
        double a = kernelSol[2 + (i1 - 1) * ncomp2 + j];
        for (p = 0; p < fwKernel * fwKernel; p++)
          effKernel[p] += a * kernel_vec[i1][p];
      }
    }

    /* Circularly shift kernel centre to (0,0) for linear convolution */
    memset(real_buf, 0, (size_t)fft_nx * fft_ny * sizeof(double));
    for (jk = 0; jk < fwKernel; jk++) {
      jj = ((jk - hwKernel) % fft_ny + fft_ny) % fft_ny;
      for (ik = 0; ik < fwKernel; ik++) {
        ii = ((ik - hwKernel) % fft_nx + fft_nx) % fft_nx;
        real_buf[ii + fft_nx * jj] = effKernel[ik + fwKernel * jk];
      }
    }

    /* FFT kernel, pointwise multiply with image spectrum, IFFT.
     * Use flat double* arithmetic: C99 _Complex and double[2] both store
     * real at offset 0 and imaginary at offset 1 (C99 §6.2.5). */
    fftw_execute_dft_r2c(plan_fwd, real_buf, ker_fft);
    {
      double* kd = (double*)ker_fft;
      const double* id = (const double*)img_fft;
      double reval, imval;
      for (p = 0; p < nc_fft; p++) {
        reval = kd[2 * p] * id[2 * p] - kd[2 * p + 1] * id[2 * p + 1];
        imval = kd[2 * p] * id[2 * p + 1] + kd[2 * p + 1] * id[2 * p];
        kd[2 * p] = reval;
        kd[2 * p + 1] = imval;
      }
    }
    fftw_execute_dft_c2r(plan_inv, ker_fft, real_buf);

    /* Normalise (FFTW does not normalise) and store valid region.
     * Extract only the valid region [0:ySize)[0:xSize) from padded real_buf
     * and store in unpadded effConv buffer. This is safe because:
     *   - real_buf is indexed with stride fft_nx (padded)
     *   - effConv[cidx] is indexed with stride xSize (unpadded)
     *   - We only access indices [0:xSize*ySize), which are allocated in both.
     * The padding region of real_buf is discarded (not needed for output).
     */
    norm = 1.0 / ((double)fft_nx * fft_ny);
    for (yp = 0; yp < ySize; yp++)
      for (xp = 0; xp < xSize; xp++)
        effConv[cidx][xp + xSize * yp] =
            (float)(real_buf[xp + fft_nx * yp] * norm);
  }

  /* --- Destroy plans and release spectral buffers --- */
#ifdef _OPENMP
#pragma omp critical(fftw_plan)
#endif
  {
    fftw_destroy_plan(plan_fwd);
    plan_fwd = NULL;
    fftw_destroy_plan(plan_inv);
    plan_inv = NULL;
  }
  fftw_free(real_buf);
  real_buf = NULL;
  fftw_free(img_fft);
  img_fft = NULL;
  fftw_free(ker_fft);
  ker_fft = NULL;
  free(effKernel);
  effKernel = NULL;

  /* --- Tile-based pixel-wise weighted sum (OpenMP parallel over tiles) ---
   * Processes image in FFT_TILE_SIZE x FFT_TILE_SIZE tiles to improve cache
   * locality and reduce L2 cache misses. For each tile, keeps all effConv[j]
   * values for pixels in that tile resident in L2 cache before proceeding to
   * next tile. Reduces prefetching overhead and improves memory bandwidth.
   * For each pixel (x,y), compute:
   *   output(x,y) = effConv[0][ni] + Σⱼ Pⱼ(xf,yf) * effConv[j+1][ni]
   * where Pⱼ = xf^ix * yf^iy in the same (ix,iy) order as make_kernel_local.
   * Each ni = xp + xSize*yp is unique; no data race on cRdata. */
  rPixX2 = 0.5 * xSize;
  rPixY2 = 0.5 * ySize;
#define FFT_TILE_SIZE 32
#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    collapse(2) private(ni, xf, yf, out, ax, ay, j, ix, iy, xp, yp, i0, j0, \
                        i2, j2) default(shared)
#endif
  for (j0 = hwKernel; j0 < ySize - hwKernel; j0 += FFT_TILE_SIZE) {
    for (i0 = hwKernel; i0 < xSize - hwKernel; i0 += FFT_TILE_SIZE) {
      int xp_max = (i0 + FFT_TILE_SIZE < xSize - hwKernel)
                       ? i0 + FFT_TILE_SIZE
                       : xSize - hwKernel;
      int yp_max = (j0 + FFT_TILE_SIZE < ySize - hwKernel)
                       ? j0 + FFT_TILE_SIZE
                       : ySize - hwKernel;

      for (yp = j0; yp < yp_max; yp++) {
        for (xp = i0; xp < xp_max; xp++) {
          ni = xp + xSize * yp;
          xf = (xp - rPixX2) / rPixX2;
          yf = (yp - rPixY2) / rPixY2;

          out = effConv[0][ni];
          j = 0;
          ax = 1.0;
          for (ix = 0; ix <= kerOrder; ix++) {
            ay = 1.0;
            for (iy = 0; iy <= kerOrder - ix; iy++) {
              out += ax * ay * effConv[j + 1][ni];
              ay *= yf;
              j++;
            }
            ax *= xf;
          }
          cRdata[ni] = (float)out;
        }
      }
    }
  }
#undef FFT_TILE_SIZE

  /* --- Free effective convolution images --- */
  for (i = 0; i < nEffConv; i++) {
    free(effConv[i]);
    effConv[i] = NULL;
  }
  free(effConv);
  effConv = NULL;

  /* --- Mask propagation + variance: block-centre direct loop ---
   * Structurally identical to the non-OMP path in spatial_convolve() with
   * the pixel accumulation (q +=) removed; cRdata is already filled above.
   * Variance is propagated exactly as in the original. */
  lkernel = (double*)calloc((size_t)fwKernel * fwKernel, sizeof(double));
  lkernel_coeffs = (double*)calloc((size_t)nCompKer, sizeof(double));
  if (!lkernel || !lkernel_coeffs) {
    LOG_ERROR("spatial_convolve_fft: out of memory (mask loop)");
    goto cleanup_fft;
  }

  nsteps_x_loc = (int)ceil((double)xSize / (double)kcStep);
  nsteps_y_loc = (int)ceil((double)ySize / (double)kcStep);

  for (j1 = 0; j1 < nsteps_y_loc; j1++) {
    j0 = j1 * kcStep + hwKernel;

    for (i = 0; i < nsteps_x_loc; i++) {
      i0 = i * kcStep + hwKernel;

      make_kernel_local_dispatch(i0 + hwKernel, j0 + hwKernel, kernelSol, lkernel,
                                 lkernel_coeffs);

      for (j2 = 0; j2 < kcStep; j2++) {
        j = j0 + j2;
        if (j >= ySize - hwKernel) break;

        for (i2 = 0; i2 < kcStep; i2++) {
          i1 = i0 + i2;
          if (i1 >= xSize - hwKernel) break;

          ni = i1 + xSize * j;
          qv = aks = uks = 0.0;
          mbit = 0x0;

          for (jc = j - hwKernel; jc <= j + hwKernel; jc++) {
            jk = j - jc + hwKernel;
            for (ic = i1 - hwKernel; ic <= i1 + hwKernel; ic++) {
              ik = i1 - ic + hwKernel;
              nc = ic + xSize * jc;
              /* Bounds check to prevent out-of-bounds access */
              if (nc < 0 || nc >= (int)image_size) {
                LOG_ERROR("spatial_convolve_fft: mask loop bounds error: nc=%d, "
                          "image_size=%zu, i1=%d, j=%d, xSize=%d",
                          nc, image_size, i1, j, xSize);
                goto cleanup_fft;
              }
              kk = lkernel[ik + jk * fwKernel];

              if (dovar) {
                if (convolveVariance)
                  qv += (*variance)[nc] * kk;
                else
                  qv += (*variance)[nc] * kk * kk;
              }

              mbit |= cMask[nc];
              aks += fabs(kk);
              if (!(cMask[nc] & FLAG_INPUT_ISBAD)) uks += fabs(kk);
            }
          }

          if (dovar) vData[ni] = (float)qv;

          mRData[ni] |= cMask[ni];
          mRData[ni] |=
              FLAG_OUTPUT_ISBAD * ((cMask[ni] & FLAG_INPUT_ISBAD) > 0);

          if (mbit) {
            if ((uks / aks) < kerFracMask)
              mRData[ni] |= (FLAG_OUTPUT_ISBAD | FLAG_BAD_CONV);
            else
              mRData[ni] |= FLAG_OK_CONV;
          }
        }
      }
    }
  }

  free(lkernel);
  lkernel = NULL;
  free(lkernel_coeffs);
  lkernel_coeffs = NULL;

  if (dovar) {
    free(*variance);
    *variance = vData;
    vData = NULL;
  }
  return;

cleanup_fft:
#ifdef _OPENMP
#pragma omp critical(fftw_plan)
#endif
{
  if (plan_fwd) {
    fftw_destroy_plan(plan_fwd);
    plan_fwd = NULL;
  }
  if (plan_inv) {
    fftw_destroy_plan(plan_inv);
    plan_inv = NULL;
  }
}
  if (real_buf) {
    fftw_free(real_buf);
    real_buf = NULL;
  }
  if (img_fft) {
    fftw_free(img_fft);
    img_fft = NULL;
  }
  if (ker_fft) {
    fftw_free(ker_fft);
    ker_fft = NULL;
  }
  if (effKernel) {
    free(effKernel);
    effKernel = NULL;
  }
  if (effConv) {
    for (i = 0; i < nEffConv; i++) free(effConv[i]);
    free(effConv);
    effConv = NULL;
  }
  if (lkernel) {
    free(lkernel);
    lkernel = NULL;
  }
  if (lkernel_coeffs) {
    free(lkernel_coeffs);
    lkernel_coeffs = NULL;
  }
  free(vData);
  vData = NULL;
}

/**
 * @brief Convolve an entire image with the spatially-varying kernel, updating
 *        the output image and propagating the pixel mask.
 *
 * @details This is the performance bottleneck of HOTPANTS (~60 % of total
 * runtime).  The image is divided into kcStep×kcStep blocks.  For each block
 * make_kernel() is called once at the block centre to evaluate the
 * spatially-varying kernel K(x,y) = Σᵢ cᵢ(x,y)·φᵢ; the same kernel image is
 * then applied to all pixels in the block by direct 2-D summation over the
 * fwKernel×fwKernel support.
 *
 * If a variance image is provided (*variance != NULL), it is propagated
 * simultaneously.  When convolveVariance is set the variance is convolved
 * linearly (appropriate for correlated noise); otherwise it is convolved with
 * K² (appropriate for independent noise).
 *
 * Mask propagation: the output mask mRData[ni] inherits the input mask of the
 * central pixel.  If any pixel within the kernel footprint is flagged and the
 * fraction of unmasked kernel weight uks/aks falls below kerFracMask, the
 * output pixel is additionally flagged FLAG_OUTPUT_ISBAD | FLAG_BAD_CONV;
 * otherwise it receives FLAG_OK_CONV.
 *
 * When built with OpenMP (-fopenmp) the outer j1 loop is parallelised.  Each
 * thread allocates its own kernel and kernel_coeffs workspace so that
 * make_kernel_local() calls are data-race free.  The mRData and cRdata writes
 * are safe because each output pixel index ni = i + xSize*j is unique across
 * iterations.
 *
 * @param image      Full-frame input image to be convolved.
 * @param variance   Pointer to the variance image pointer; updated in-place if
 *                   non-NULL (old allocation freed and replaced).
 * @param xSize      Image width in pixels.
 * @param ySize      Image height in pixels.
 * @param kernelSol  Kernel solution vector produced by fitKernel().
 * @param cRdata     Output: full-frame convolved image (pre-allocated by
 *                   caller).
 * @param cMask      Input pixel mask array used for mask propagation.
 */
/**
 * @brief Apply spatially-varying convolution kernel to an image region.
 *
 * @details Evaluates K(x,y) = Σᵢ cᵢ(x,y)·φᵢ at each pixel and computes the
 * convolution D(x,y) = I(x,y) - [T ⊗ K](x,y). Two implementations available:
 * - Direct convolution: O(N·k²) pixel-level loops (fallback)
 * - FFT acceleration: O(N log N) precomputed basis convolutions (preferred for
 * large images)
 *
 * Dispatches to spatial_convolve_fft() if USE_FFTW is enabled and FFTW3 is
 * available; otherwise uses direct O(N·k²) implementation.
 *
 * @param[in] image      Flat image array (e.g., template) to be convolved.
 * @param[in,out] variance Optional variance map (pointer-to-pointer); if
 * non-NULL, noise propagation is computed in-place.
 * @param[in] xSize      Image width in pixels.
 * @param[in] ySize      Image height in pixels.
 * @param[in] kernelSol  Kernel solution vector cᵢ from fitKernel();
 * dimensionality nCompKer.
 * @param[out] cRdata    Output convolved image (pre-allocated by caller;
 * xSize×ySize).
 * @param[in] cMask      Input pixel mask for mask propagation (0 = bad,
 * non-zero = good).
 *
 * @note Calls make_kernel() to evaluate K(x,y) at each pixel; expensive
 * operation. For multi-threaded execution, loop indices and temporary buffers
 * are thread-local.
 * @see spatial_convolve_fft() for FFT-accelerated variant
 * @see fitKernel() for kernel solution computation
 */
void spatial_convolve(float* image, float** variance, int xSize, int ySize,
                      double* kernelSol, float* cRdata, int* cMask) {
  spatial_convolve_fft(image, variance, xSize, ySize, kernelSol, cRdata,
                       cMask);
}

/**
 * @brief Dispatcher: evaluate kernel at (xi, yi) using polynomial or TPS.
 *
 * @details Routes to make_kernel() or make_kernel_tps() based on useTPS flag.
 * This is the primary entry point for kernel evaluation during convolution.
 *
 * @param xi         x pixel coordinate
 * @param yi         y pixel coordinate
 * @param kernelSol  Kernel solution vector
 * @return Kernel sum (photometric scaling factor)
 *
 * @see make_kernel() for polynomial evaluation
 * @see make_kernel_tps() for TPS evaluation
 */
double make_kernel_dispatch(int xi, int yi, double* kernelSol) {
  /* Dispatch based on kernel basis type (Gaussian vs. Delta) and spatial variation mode
     (polynomial vs. TPS).
     NOTE: Phase 2+ will add iBasisType == BASIS_TYPE_DELTA branch. */
  if (iBasisType == BASIS_TYPE_DELTA) {
    LOG_ERROR("Delta function basis not yet implemented (Phase 2+)");
    return 0.0;
    /* TODO: return make_kernel_delta(xi, yi, kernelSol); */
  }

  /* Gaussian basis (default) */
  if (useTPS) {
    return make_kernel_tps(xi, yi, kernelSol);
  } else {
    return make_kernel(xi, yi, kernelSol);
  }
}

/**
 * @brief Evaluate the spatially-varying convolution kernel at image position
 *        (xi, yi) and return its pixel sum.
 *
 * @details Computes the spatial polynomial weights at (xi, yi) (normalised to
 * [-1, 1] across the image) and uses them to form the linear combination of
 * kernel coefficients:
 *   kernel_coeffs[i1] = Σ_{ix,iy} kernelSol[k] · xf^ix · yf^iy
 * then assembles the kernel image:
 *   kernel[p] = Σ_{i1} kernel_coeffs[i1] · kernel_vec[i1][p]
 * The result is stored in the global array kernel[] and its pixel sum is
 * returned.  The kernel sum equals the photometric flux-scaling ratio between
 * convolved and reference images at this position.
 *
 * @param xi         x pixel coordinate at which to evaluate the kernel.
 * @param yi         y pixel coordinate at which to evaluate the kernel.
 * @param kernelSol  Kernel solution vector (output of fitKernel()).
 * @return Sum of all pixels in the assembled kernel image (the flux-scaling
 *         factor at this position).
 */
double make_kernel(int xi, int yi, double* kernelSol) {
  int gaussianCompIdx, solutionIdx, polyDegX, polyDegY, kernelPixelIdx;
  int pixelCompIdx;
  double polyBasisX, polyBasisY, kernelSum;
  double normalizedX, normalizedY;

  solutionIdx = 2;
  /* RANGE FROM -1 to 1 */
  normalizedX = (xi - 0.5 * rPixX) / (0.5 * rPixX);
  normalizedY = (yi - 0.5 * rPixY) / (0.5 * rPixY);

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
 * @brief Evaluate the spatially-varying convolution kernel using thin plate
 * spline interpolation.
 *
 * @details TPS-based alternative to make_kernel(). For each kernel component,
 * evaluates the fitted TPS RBF surface at (xi, yi) to obtain the spatially-
 * varying coefficient. Then assembles the final kernel image as a weighted
 * sum of basis functions.
 *
 * TPS parameters are stored in kernelSol when useTPS=1. Layout:
 *   [0..poly_size-1]: polynomial coefficients (for background compat)
 *   [poly_size..poly_size+nCompKer*nS-1]: RBF weights
 *   [poly_size+nCompKer*nS..poly_size+nCompKer*nS+3*nCompKer-1]: poly trends
 *   [poly_size+nCompKer*nS+3*nCompKer..end]: stamp positions
 *
 * @param xi         x pixel coordinate at which to evaluate
 * @param yi         y pixel coordinate at which to evaluate
 * @param kernelSol  Extended kernel solution vector (includes TPS parameters)
 * @return Sum of all pixels in the assembled kernel image
 *
 * @see tps_evaluate() for RBF surface evaluation
 * @see make_kernel() for polynomial (non-TPS) version
 */
double make_kernel_tps(int xi, int yi, double* kernelSol) {
  int gaussianCompIdx, kernelPixelIdx, pixelCompIdx;
  double kernelSum;
  double *tps_weights, *tps_poly, *positions;
  int pos_offset;

  /* For TPS evaluation, nS is the number of stamps per region used during fit.
     Note: This assumes single-region (nRegX=1, nRegY=1), so nS == nStampX*nStampY
  */

  pos_offset = kernelSol_offset_tps_positions(nS);

  /* Evaluate TPS-interpolated kernel coefficients at (xi, yi) */
  for (gaussianCompIdx = 1; gaussianCompIdx < nCompKer; gaussianCompIdx++) {
    tps_weights = kernelSol + kernelSol_offset_tps_weights(gaussianCompIdx, nS);
    tps_poly = kernelSol + kernelSol_offset_tps_poly(gaussianCompIdx, nS);
    positions = kernelSol + pos_offset;

    kernel_coeffs[gaussianCompIdx] =
        tps_evaluate((double)xi, (double)yi, positions, nS, tps_weights, tps_poly);
  }

  /* Constant component (no spatial variation) */
  kernel_coeffs[0] = kernelSol[1];

  /* Assemble kernel image from basis functions */
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
 * @brief Evaluate the spatially-varying background using thin plate spline
 * interpolation.
 *
 * @details TPS-based alternative to polynomial get_background(). Evaluates the
 * fitted TPS RBF surface at (xi, yi) to obtain the background value.
 *
 * @param xi         x pixel coordinate.
 * @param yi         y pixel coordinate.
 * @param kernelSol  Extended kernel solution vector (includes background TPS parameters)
 * @return Background value at (xi, yi)
 *
 * @see get_background() for polynomial (non-TPS) version
 */
double get_background_tps(int xi, int yi, double* kernelSol) {
  double *tps_weights, *tps_poly, *positions;
  int pos_offset, n_stamps;

  /* nS is the number of stamps per region used during fit */
  n_stamps = nS;
  pos_offset = kernelSol_offset_tps_positions(n_stamps);

  tps_weights = kernelSol + kernelSol_offset_tps_bg_weights(n_stamps);
  tps_poly = kernelSol + kernelSol_offset_tps_bg_poly(n_stamps);
  positions = kernelSol + pos_offset;

  return tps_evaluate((double)xi, (double)yi, positions, n_stamps, tps_weights, tps_poly);
}

/**
 * @brief Evaluate the spatially-varying background polynomial at image position
 *        (xi, yi).
 *
 * @details The background model is a 2-D polynomial of order bgOrder, with
 * coefficients stored in kernelSol starting at index ncompBG+1.  Coordinates
 * are normalised to [-1, 1].  The returned value is the additive sky background
 * difference between convolved and reference images at (xi, yi), to be
 * subtracted when forming the difference image.
 *
 * @param xi         x pixel coordinate.
 * @param yi         y pixel coordinate.
 * @param kernelSol  Kernel solution vector; background coefficients begin
 *                   immediately after the kernel coefficient block.
 * @return Background value at (xi, yi) in data units.
 */
double get_background(int xi, int yi, double* kernelSol) {
  if (useTPS) {
    return get_background_tps(xi, yi, kernelSol);
  }

  double background, polyBasisX, polyBasisY, normalizedX, normalizedY;
  int polyDegX, polyDegY, solutionIdx;
  int backgroundComponentOffset;

  backgroundComponentOffset =
      (nCompKer - 1) * (((kerOrder + 1) * (kerOrder + 2)) / 2) + 1;

  background = 0.0;
  solutionIdx = 1;
  /* RANGE FROM -1 to 1 */
  normalizedX = (xi - 0.5 * rPixX) / (0.5 * rPixX);
  normalizedY = (yi - 0.5 * rPixY) / (0.5 * rPixY);

  polyBasisX = 1.0;
  for (polyDegX = 0; polyDegX <= bgOrder; polyDegX++) {
    polyBasisY = 1.0;
    for (polyDegY = 0; polyDegY <= bgOrder - polyDegX; polyDegY++) {
      background += kernelSol[backgroundComponentOffset + solutionIdx++] *
                    polyBasisX * polyBasisY;
      /* fprintf(stderr, "bg: %d %d %d %d %f %f %f\n", xi, yi, polyDegX,
       * polyDegY, polyBasisX, polyBasisY,
       * kernelSol[backgroundComponentOffset+solutionIdx-1]); */
      polyBasisY *= normalizedY;
    }
    polyBasisX *= normalizedX;
  }
  return background;
}

/**
 * @brief Construct the model convolved image for a substamp from the kernel
 *        solution and the precomputed convolution vectors.
 *
 * @details Evaluates the spatial polynomial weights at the substamp centre and
 * forms the model prediction:
 *   csModel[p] = Σ_{i1} coeff[i1] · stamp->vectors[i1][p]
 * where coeff[i1] is the spatially-evaluated kernel coefficient for basis i1.
 * The result is the predicted appearance of the template convolved to match the
 * reference at this substamp location, and is used by getStampSig() to compute
 * residual metrics.  Background is not included here; it is added separately
 * via get_background().
 *
 * @param stamp      Stamp with pre-filled convolution vectors.
 * @param kernelSol  Kernel solution vector.
 * @param csModel    Output: fwKSStamp×fwKSStamp float array filled with the
 *                   model prediction (pre-allocated by caller).
 */
void make_model(stamp_struct* stamp, double* kernelSol, float* csModel) {
  int gaussianCompIdx, solutionIdx, polyDegX, polyDegY, stampPixelIdx,
      stampCenterX, stampCenterY;
  double polyBasisX, polyBasisY, polynomialCoeff;
  double* vector;
  double normalizedX, normalizedY;

  stampCenterX = stamp->xss[stamp->sscnt];
  stampCenterY = stamp->yss[stamp->sscnt];

  /* RANGE FROM -1 to 1 */
  normalizedX = (stampCenterX - 0.5 * rPixX) / (0.5 * rPixX);
  normalizedY = (stampCenterY - 0.5 * rPixY) / (0.5 * rPixY);

  for (stampPixelIdx = 0; stampPixelIdx < fwKSStamp * fwKSStamp;
       stampPixelIdx++)
    csModel[stampPixelIdx] = 0.0;

  vector = stamp->vectors[0];
  polynomialCoeff = kernelSol[1];
  for (stampPixelIdx = 0; stampPixelIdx < fwKSStamp * fwKSStamp;
       stampPixelIdx++)
    csModel[stampPixelIdx] += polynomialCoeff * vector[stampPixelIdx];

  solutionIdx = 2;
  for (gaussianCompIdx = 1; gaussianCompIdx < nCompKer; gaussianCompIdx++) {
    vector = stamp->vectors[gaussianCompIdx];
    polynomialCoeff = 0.0;
    polyBasisX = 1.0;
    for (polyDegX = 0; polyDegX <= kerOrder; polyDegX++) {
      polyBasisY = 1.0;
      for (polyDegY = 0; polyDegY <= kerOrder - polyDegX; polyDegY++) {
        polynomialCoeff += kernelSol[solutionIdx++] * polyBasisX * polyBasisY;
        polyBasisY *= normalizedY;
      }
      polyBasisX *= normalizedX;
    }

    for (stampPixelIdx = 0; stampPixelIdx < fwKSStamp * fwKSStamp;
         stampPixelIdx++) {
      csModel[stampPixelIdx] += polynomialCoeff * vector[stampPixelIdx];
    }
  }
  return;
}
