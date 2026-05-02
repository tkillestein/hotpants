#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include<fitsio.h>
#include<lapacke.h>
#ifdef USE_FFTW
#include<fftw3.h>
#endif

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

/**
 * @brief Precompute all kernel basis function images and store them in kernel_vec.
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
    int gaussIdx, idegx, idegy, nvec;
    int ren;

    nvec = 0;
    for (gaussIdx = 0; gaussIdx < ngauss; gaussIdx++) {
        for (idegx = 0; idegx <= deg_fixe[gaussIdx]; idegx++) {
            for (idegy = 0; idegy <= deg_fixe[gaussIdx]-idegx; idegy++) {
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
int fillStamp(stamp_struct *stamp, float *imConv, float *imRef) {

    int       renormalizeFlag = 0;
    int       pixelX,pixelY,stampCenterX,stampCenterY,pixelOffsetX,pixelOffsetY,idegx,idegy,bgDegX,bgDegY,substampIndexX,substampIndexY,vectorCompIdx,gaussianCompIdx,vectorComponentIdx;
    double    polyBasisX,polyBasisY,normalizedX,normalizedY;
    double *im;
    float     halfPixX, halfPixY;

    halfPixX   = 0.5 * rPixX;
    halfPixY   = 0.5 * rPixY;

    if (verbose >= 1)
        fprintf(stderr, "    xs  : %4i ys  : %4i sig: %6.3f sscnt: %4i nss: %4i \n",
                stamp->x, stamp->y, stamp->chi2, stamp->sscnt, stamp->nss);
    if (stamp->sscnt >= stamp->nss) {
        /* have gone through all the good substamps, reject this stamp */
        /*if (verbose >= 2) fprintf(stderr, "    ******** REJECT stamp (out of substamps)\n");*/
        if (verbose >= 1)
            fprintf(stderr, "        Reject stamp\n");
        return 1;
    }

    /* ====================================================================
       POLYNOMIAL BASIS CONSTRUCTION
       ====================================================================
       Iterate over kernel basis triplets (gaussianCompIdx, idegx, idegy):
       - gaussianCompIdx: Gaussian component index [0, ngauss)
       - idegx, idegy: polynomial degrees such that idegx + idegy <= deg_fixe[gaussianCompIdx]

       This implements the triangular basis expansion:
         K(x,y) = Sum_i c_i(x,y) * phi_i
       where phi_i = exp(-(x^2+y^2)*sigma^2) * x^deg_x * y^deg_y

       The counter vectorComponentIdx sequentially assigns indices to each (gaussianCompIdx, idegx, idegy)
       triplet; these indices are used by xy_conv_stamp() to access precomputed
       convolution responses stored in stamp->vectors[vectorComponentIdx].

       Alard & Lupton (1998), Sect. 2: "The kernel is expanded as a sum of
       Gaussian basis elements with polynomial spatial weighting."
       Reference: https://iopscience.iop.org/article/10.1086/305984
       ==================================================================== */
    vectorComponentIdx = 0;
    for (gaussianCompIdx = 0; gaussianCompIdx < ngauss; gaussianCompIdx++) {
        for (idegx = 0; idegx <= deg_fixe[gaussianCompIdx]; idegx++) {
            for (idegy = 0; idegy <= deg_fixe[gaussianCompIdx]-idegx; idegy++) {

                renormalizeFlag = 0;
                /* Detect odd-degree terms: (deg/2)*2 - deg evaluates to
                   -1 if deg is odd, 0 if even. Used to trigger renormalization
                   of even-parity basis functions (those with pixelOffsetX==0 && pixelOffsetY==0). */
                pixelOffsetX = (idegx / 2) * 2 - idegx;
                pixelOffsetY = (idegy / 2) * 2 - idegy;
                if (pixelOffsetX == 0 && pixelOffsetY == 0 && vectorComponentIdx > 0)
                    renormalizeFlag = 1;  /* Set renormalization flag for higher-order even-parity basis */

                /* fill stamp->vectors[vectorComponentIdx] with convolved image */
                /* image is convolved with functional form of kernel, fit later for amplitude */
                xy_conv_stamp(stamp, imConv, vectorComponentIdx, renormalizeFlag);
                ++vectorComponentIdx;
            }
        }
    }

    /* get the krefArea data */
    if (cutSStamp(stamp, imRef))
        return 1;

    /* fill stamp->vectors[vectorComponentIdx+++] with x^(bg) * y^(bg) for background fit */
    stampCenterX = stamp->xss[stamp->sscnt];
    stampCenterY = stamp->yss[stamp->sscnt];
    substampIndexX = stampCenterX - hwKSStamp;
    substampIndexY = stampCenterY - hwKSStamp;
    for (pixelX = stampCenterX - hwKSStamp; pixelX <= stampCenterX + hwKSStamp; pixelX++) {
        normalizedX = (pixelX - halfPixX) / halfPixX;

        for (pixelY = stampCenterY - hwKSStamp; pixelY <= stampCenterY + hwKSStamp; pixelY++) {
            /* fprintf(stderr, "%d %d %d %d %d %d\n", k, stampCenterX, stampCenterY,pixelX, pixelY, fwKSStamp); */
            normalizedY = (pixelY - halfPixY) / halfPixY;

            polyBasisX = 1.0;
            vectorCompIdx = vectorComponentIdx;
            for (bgDegX = 0; bgDegX <= bgOrder; bgDegX++) {
                polyBasisY = 1.0;
                for (bgDegY = 0; bgDegY <= bgOrder - bgDegX; bgDegY++) {
                    im = stamp->vectors[vectorCompIdx];
                    im[pixelX-substampIndexX+fwKSStamp*(pixelY-substampIndexY)] = polyBasisX * polyBasisY;
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
 * If usePCA is set the work is delegated to kernel_vector_PCA().
 *
 * @param n     Sequential index of this basis function within kernel_vec.
 * @param deg_x Polynomial degree of the x-direction filter.
 * @param deg_y Polynomial degree of the y-direction filter.
 * @param ig    Gaussian component index (selects sigma_gauss[ig]).
 * @param ren   Output flag: set to 1 when the zeroth basis is subtracted.
 * @return Pointer to a newly allocated fwKernel*fwKernel double array; the
 *         caller stores this in kernel_vec[n] and must eventually free it.
 */
double *kernel_vector(int n, int deg_x, int deg_y, int ig, int *ren) {
    double    *vector=NULL,*kernel0=NULL;
    int       kernelRow,kernelCol,xFilterIdx,xParityCheck,yParityCheck,xKernelIdx,pixelIdx;
    double    filterXSum,filterYSum,kernelPosition,gaussianVal;

    if (usePCA) {
        return kernel_vector_PCA(n, deg_x, deg_y, ig, ren);
    }

    vector = (double *)malloc(fwKernel*fwKernel*sizeof(double));
    /* Detect parity: (deg/2)*2 - deg = 0 if even, -1 if odd.
       Odd-degree filters (xParityCheck != 0 || yParityCheck != 0) retain full normalization;
       even-degree filters (xParityCheck==0 && yParityCheck==0) are normalized and (for n>0)
       have the n=0 basis subtracted to form differential corrections. */
    xParityCheck = (deg_x / 2) * 2 - deg_x;
    yParityCheck = (deg_y / 2) * 2 - deg_y;
    filterXSum = filterYSum = 0.0;
    *ren = 0;

    for (xKernelIdx = 0; xKernelIdx < fwKernel; xKernelIdx++) {
        kernelPosition            = (double)(xKernelIdx - hwKernel);
        xFilterIdx            = xKernelIdx+n*fwKernel;
        gaussianVal           = exp(-kernelPosition * kernelPosition * sigma_gauss[ig]);
        filter_x[xFilterIdx]  = gaussianVal * pow(kernelPosition, deg_x);
        filter_y[xFilterIdx]  = gaussianVal * pow(kernelPosition, deg_y);
        filterXSum       += filter_x[xFilterIdx];
        filterYSum       += filter_y[xFilterIdx];
    }

    if (n > 0)
        kernel0 = kernel_vec[0];

    filterXSum = 1. / filterXSum;
    filterYSum = 1. / filterYSum;

    if (xParityCheck == 0 && yParityCheck == 0) {
        for (xKernelIdx = 0; xKernelIdx < fwKernel; xKernelIdx++) {
            filter_x[xKernelIdx+n*fwKernel] *= filterXSum;
            filter_y[xKernelIdx+n*fwKernel] *= filterYSum;
        }

        for (kernelRow = 0; kernelRow < fwKernel; kernelRow++) {
            for (kernelCol = 0; kernelCol < fwKernel; kernelCol++) {
                vector[kernelRow+fwKernel*kernelCol] = filter_x[kernelRow+n*fwKernel] * filter_y[kernelCol+n*fwKernel];
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
                vector[kernelRow+fwKernel*kernelCol] = filter_x[kernelRow+n*fwKernel] * filter_y[kernelCol+n*fwKernel];
            }
        }
    }
    return vector;
}

/**
 * @brief Compute the n-th kernel basis function image from a pre-supplied PCA
 *        basis set.
 *
 * @details Instead of Gaussian-polynomial analytic forms, each basis image is
 * read directly from the global PCA[n] array (empirical principal components
 * derived from observed PSFs).  As with kernel_vector(), for n > 0 the zeroth
 * PCA component is subtracted so that higher components represent differential
 * corrections, and *ren is set to 1 to propagate this renormalisation flag to
 * xy_conv_stamp_PCA().
 *
 * @param n     Sequential index of this basis function.
 * @param deg_x Polynomial degree (unused in PCA mode; kept for API symmetry).
 * @param deg_y Polynomial degree (unused in PCA mode; kept for API symmetry).
 * @param ig    Gaussian index (unused in PCA mode; kept for API symmetry).
 * @param ren   Output flag: set to 1 when the zeroth PCA basis is subtracted.
 * @return Pointer to a newly allocated fwKernel*fwKernel double array
 *         initialised with PCA[n] (minus PCA[0] for n > 0).
 */
double *kernel_vector_PCA(int n, int deg_x, int deg_y, int ig, int *ren) {
    double    *vector=NULL,*kernel0=NULL;
    int xIdx,yIdx;

    vector = (double *)malloc(fwKernel*fwKernel*sizeof(double));

    for (xIdx = 0; xIdx < fwKernel; xIdx++) {
        for (yIdx = 0; yIdx < fwKernel; yIdx++) {
            vector[xIdx+fwKernel*yIdx] = PCA[n][xIdx+fwKernel*yIdx];
        }
    }

    if (n > 0)
        kernel0 = kernel_vec[0];

    if (n > 0) {
        for (xIdx = 0; xIdx < fwKernel * fwKernel; xIdx++) {
            vector[xIdx] -= kernel0[xIdx];
        }
        *ren = 1;
    }

    return vector;
}

/**
 * @brief Convolve the substamp region of an image with the n-th separable
 *        Gaussian-polynomial kernel basis function.
 *
 * @details **Separability optimization:** Gaussian-polynomial filters are separable:
 *   φ(x,y) = filter_x[x] * filter_y[y]
 * where filter_x[k] = exp(-k^2*σ^2) * k^deg_x and similarly for filter_y.
 *
 * Standard 2D convolution costs O(k^2*n^2) FLOPs (k = kernel width, n = image dimension).
 * Separable 1D passes reduce this to O(2*k*n^2), a speedup factor of k/2 (typically
 * 5–10×). This is the primary reason HOTPANTS achieves real-time performance on
 * large images.
 *
 * The implementation uses the separability of the Gaussian-polynomial filter to perform
 * the convolution in two 1-D passes (first along y, then along x), storing the
 * fwKSStamp×fwKSStamp result in stamp->vectors[n].  This is the performance-
 * critical inner loop (~60 % of total runtime according to profiling).
 *
 * If the renormalise flag @p ren is set (i.e. kernel_vector() returned ren=1
 * for this basis), the zeroth-basis convolved image stamp->vectors[0] is
 * subtracted pixel-by-pixel from the result, matching the differential
 * construction of the higher-order basis images.
 *
 * In PCA mode the call is forwarded to xy_conv_stamp_PCA().
 *
 * Reference: Alard & Lupton (1998), Sect. 2; classic signal-processing optimization.
 *
 * @param stamp  Stamp whose current substamp (stamp->sscnt) defines the pixel
 *               region to convolve.
 * @param image  Full-frame input image to be convolved.
 * @param n      Index of the kernel basis function; selects filter_x/filter_y
 *               rows and the output slot stamp->vectors[n].
 * @param ren    If non-zero, subtract stamp->vectors[0] from the result after
 *               convolution (renormalisation step for higher-order bases).
 */
void xy_conv_stamp(stamp_struct *stamp, float *image, int n, int ren) {
    int       stampColIdx,stampRowIdx,kerOffsetX,kerOffsetY,xij,sub_width,xi,yi,imgColIdx,imgRowIdx,pixelIdx;
    double    *v0,*imc;


    if (usePCA) {
        xy_conv_stamp_PCA(stamp, image, n, ren);
        return;
    }

    xi  = stamp->xss[stamp->sscnt];
    yi  = stamp->yss[stamp->sscnt];
    imc = stamp->vectors[n];

    sub_width = fwKSStamp + fwKernel - 1;

    /* pull area to convolve out of full reference image region */
    /* convolve with y filter */
    for(imgColIdx = xi - hwKSStamp - hwKernel; imgColIdx <= xi + hwKSStamp + hwKernel; imgColIdx++) {
        for(imgRowIdx = yi - hwKSStamp; imgRowIdx <= yi + hwKSStamp; imgRowIdx++) {
            xij = imgColIdx - xi + sub_width / 2 + sub_width * (imgRowIdx - yi + hwKSStamp);
            temp[xij] = 0.0;
            for(kerOffsetY = -hwKernel; kerOffsetY <= hwKernel; kerOffsetY++) {
                temp[xij] += image[imgColIdx+rPixX*(imgRowIdx+kerOffsetY)] * filter_y[hwKernel-kerOffsetY+n*fwKernel];
            }
        }
    }

    /* convolve with x filter */
    for(stampRowIdx = -hwKSStamp; stampRowIdx <= hwKSStamp; stampRowIdx++) {
        for(stampColIdx = -hwKSStamp; stampColIdx <= hwKSStamp;stampColIdx++) {
            xij = stampColIdx + hwKSStamp + fwKSStamp * (stampRowIdx + hwKSStamp);
            imc[xij] = 0.0;
            for(kerOffsetX = -hwKernel; kerOffsetX <= hwKernel; kerOffsetX++) {
                imc[xij] += temp[stampColIdx+kerOffsetX+sub_width/2+sub_width*(stampRowIdx+hwKSStamp)] * filter_x[hwKernel-kerOffsetX+n*fwKernel];
            }
        }
    }

    if (ren) {
        v0 = stamp->vectors[0];
        for(pixelIdx = 0; pixelIdx < fwKSStamp * fwKSStamp; pixelIdx++) imc[pixelIdx] -= v0[pixelIdx];
    }
    
    return;
}

/**
 * @brief Convolve the substamp region of an image with the n-th PCA kernel
 *        basis function using a direct 2-D summation.
 *
 * @details Unlike xy_conv_stamp(), which exploits filter separability, this
 * routine performs a full 2-D correlation between the image patch and the
 * fwKernel×fwKernel PCA basis image PCA[n].  The fwKSStamp×fwKSStamp result
 * is written into stamp->vectors[n].  If @p ren is non-zero, stamp->vectors[0]
 * is subtracted to match the differential construction used in PCA mode.
 *
 * @param stamp  Stamp whose current substamp defines the pixel region to convolve.
 * @param image  Full-frame input image to be convolved.
 * @param n      PCA basis index; selects PCA[n] and the output slot
 *               stamp->vectors[n].
 * @param ren    If non-zero, subtract stamp->vectors[0] from the result after
 *               convolution.
 */
void xy_conv_stamp_PCA(stamp_struct *stamp, float *image, int n, int ren) {

    int       imgColIdx,imgRowIdx,kerOffsetX,kerOffsetY,xij,xi,yi,pixelIdx;
    double    *v0,*imc;

    xi  = stamp->xss[stamp->sscnt];
    yi  = stamp->yss[stamp->sscnt];
    imc = stamp->vectors[n];

    /* pull area to convolve out of full reference image region */
    for(imgRowIdx = yi - hwKSStamp; imgRowIdx <= yi + hwKSStamp; imgRowIdx++) {
        for(imgColIdx = xi - hwKSStamp; imgColIdx <= xi + hwKSStamp; imgColIdx++) {
            xij      = imgColIdx - (xi - hwKSStamp) + fwKSStamp * (imgRowIdx - (yi - hwKSStamp));
            imc[xij] = 0.;

            for(kerOffsetY = -hwKernel; kerOffsetY <= hwKernel; kerOffsetY++) {
                for(kerOffsetX = -hwKernel; kerOffsetX <= hwKernel; kerOffsetX++) {
                    imc[xij] += image[(imgColIdx+kerOffsetX)+rPixX*(imgRowIdx+kerOffsetY)] * PCA[n][(kerOffsetX+hwKernel) + fwKernel*(kerOffsetY+hwKernel)];
                }
            }
        }
    }

    if (ren) {
        v0 = stamp->vectors[0];
        for(pixelIdx = 0; pixelIdx < fwKSStamp * fwKSStamp; pixelIdx++) imc[pixelIdx] -= v0[pixelIdx];
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
static int solve_spd(double **a, int n, double *b) {
    lapack_int info;
    int i, j;
    double *flat = (double *)malloc((size_t)n * n * sizeof(double));
    double *rhs  = (double *)malloc((size_t)n * sizeof(double));
    if (!flat || !rhs) { free(flat); free(rhs); return 1; }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            flat[i * n + j] = a[i + 1][j + 1];
        rhs[i] = b[i + 1];
    }

    info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', (lapack_int)n, flat, (lapack_int)n);
    if (info == 0)
        info = LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'U', (lapack_int)n, 1,
                              flat, (lapack_int)n, rhs, 1);
    if (info == 0)
        for (i = 0; i < n; i++)
            b[i + 1] = rhs[i];
    else
        fprintf(stderr, "WARNING: solve_spd failed (LAPACKE info=%d); "
                "kernel solution may be unreliable\n", (int)info);

    free(flat);
    free(rhs);
    return (info != 0) ? 1 : 0;
}

/**
 * @brief Solve for the spatially-varying kernel by fitting the global normal
 *        equations, iterating until stamp sigma-clipping converges.
 *
 * @details Implements the full Alard & Lupton (1998) kernel solution:
 *  1. Assembles the global normal-equations matrix M via build_matrix() and
 *     the right-hand-side vector b via build_scprod().
 *  2. Solves M·a = b by LU decomposition (ludcmp / lubksb); the solution
 *     vector a is written into @p kernelSol.
 *  3. Calls check_again() to evaluate each stamp's residual sigma and reject
 *     or replace bad substamps; if any stamp is flagged the matrices are
 *     rebuilt and the linear system is solved again.
 *  4. Iterates steps 1–3 until check_again() reports no further changes.
 *
 * Note: the normal-equations matrix is symmetric positive-definite; the
 * correct solver is Cholesky decomposition, but LU is used here.
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
void fitKernel(stamp_struct *stamps, float *imRef, float *imConv, float *imNoise, double *kernelSol,
               double *meansigSubstamps, double *scatterSubstamps, int *NskippedSubstamps) {

    double **matrix;
    char convergedFlag;
    int matrixRowIdx,mat_size;
    int ncomp1, ncomp2, ncomp, nbg_vec;

    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;

    mat_size   = ncomp1 * ncomp2 + nbg_vec + 1;

    if (verbose >= 2) fprintf(stderr, " Mat_size: %i ncomp2: %i ncomp1: %i nbg_vec: %i \n",
                              mat_size, ncomp2, ncomp1, nbg_vec);

    /* allocate fitting matrix */
    matrix = (double **)malloc((mat_size + 1)*sizeof(double *));
    for (matrixRowIdx = 0; matrixRowIdx <= mat_size; matrixRowIdx++)
        matrix[matrixRowIdx] = (double *)malloc((mat_size + 1)*sizeof(double));

    /* allocate weight matrix */
    wxy = (double **)malloc(nS*sizeof(double *));
    for (matrixRowIdx = 0; matrixRowIdx < nS; matrixRowIdx++)
        wxy[matrixRowIdx] = (double *)malloc(ncomp2*sizeof(double));



    if (verbose>=2) fprintf(stderr, " Expanding Matrix For Full Fit\n");
    build_matrix(stamps, nS, matrix);
    build_scprod(stamps, nS, imRef, kernelSol);

    solve_spd(matrix, mat_size, kernelSol);

    if (verbose>=2) fprintf(stderr, " Checking again\n");
    convergedFlag = check_again(stamps, kernelSol, imConv, imRef, imNoise, meansigSubstamps, scatterSubstamps, NskippedSubstamps);

    while(convergedFlag) {

        fprintf(stderr, "\n Re-Expanding Matrix\n");
        build_matrix(stamps, nS, matrix);
        build_scprod(stamps, nS, imRef, kernelSol);

        solve_spd(matrix, mat_size, kernelSol);

        fprintf(stderr, " Checking again\n");
        convergedFlag = check_again(stamps, kernelSol, imConv, imRef, imNoise, meansigSubstamps, scatterSubstamps, NskippedSubstamps);
    }
    fprintf(stderr, " Sigma clipping of bad stamps converged, kernel determined\n");

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
void build_matrix0(stamp_struct *stamp) {

    int       kernelBasisIdx1,kernelBasisIdx2,pixStamp,pixelIdx,kernelBasisIdx,bgVectorIdx=0;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    p0,q;
    double    **vec;

    ncomp1   = nCompKer;
    ncomp2   = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp    = ncomp1 * ncomp2;
    nbg_vec  = ((bgOrder + 1) * (bgOrder + 2)) / 2;

    pixStamp = fwKSStamp * fwKSStamp;

    vec      = stamp->vectors;

    /* loop over the convolved images created by xy_conv_stamp() */
    /* each level represents ngauss and deg_gauss */
    for (kernelBasisIdx1 = 0; kernelBasisIdx1 < ncomp1; kernelBasisIdx1++) {
        for (kernelBasisIdx2 = 0; kernelBasisIdx2 <= kernelBasisIdx1; kernelBasisIdx2++) {
            q = 0.0;
            /* integrate W_m1 and W_m2 (sum over all pixels) */
            for (pixelIdx = 0; pixelIdx < pixStamp; pixelIdx++)
                q += vec[kernelBasisIdx1][pixelIdx] * vec[kernelBasisIdx2][pixelIdx];

            /* Q from Eqn 3. in Alard */
            stamp->mat[kernelBasisIdx1+1][kernelBasisIdx2+1] = q;
        }
    }

    for (kernelBasisIdx = 0; kernelBasisIdx < ncomp1; kernelBasisIdx++) {
        /* bgVectorIdx = index into the background vector array (the constant background term) */
        bgVectorIdx = ncomp1;

        p0 = 0.0;
        /* integrate convolved images and first order background (equals 1 everywhere!)*/
        for (pixelIdx = 0; pixelIdx < pixStamp; pixelIdx++)
            p0 += vec[kernelBasisIdx][pixelIdx] * vec[bgVectorIdx][pixelIdx];
        stamp->mat[ncomp1+1][kernelBasisIdx+1] = p0;
    }

    /* integrate first order background with itself */
    /* NOTE : DON'T MASK K HERE - BACKGROUND! */
    for (pixelIdx = 0, q = 0.0; pixelIdx < pixStamp; pixelIdx++)
        q += vec[bgVectorIdx][pixelIdx] * vec[ncomp1][pixelIdx];
    stamp->mat[ncomp1+1][ncomp1+1] = q;

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
void build_scprod0(stamp_struct *stamp, float *image) {

    int       stampColIdx,stampRowIdx,xi,yi,kernelBasisIdx,pixelIdx;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    p0,q;
    double **vec;

    ncomp1  = nCompKer;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;

    vec = stamp->vectors;
    xi  = stamp->xss[stamp->sscnt];
    yi  = stamp->yss[stamp->sscnt];

    /* Do eqn 4. in Alard */

    /* Multiply each order's convolved image with reference image */
    for (kernelBasisIdx = 0; kernelBasisIdx < ncomp1; kernelBasisIdx++) {
        p0 = 0.0;
        for (stampColIdx = -hwKSStamp; stampColIdx <= hwKSStamp; stampColIdx++) {
            for (stampRowIdx = -hwKSStamp; stampRowIdx <= hwKSStamp; stampRowIdx++) {
                pixelIdx   = stampColIdx + hwKSStamp + fwKSStamp * (stampRowIdx + hwKSStamp);
                p0 += vec[kernelBasisIdx][pixelIdx] * image[stampColIdx+xi+rPixX*(stampRowIdx+yi)];
            }
        }
        stamp->scprod[kernelBasisIdx+1] = p0;
    }

    /* Multiply first order background model with reference image */
    q = 0.0;
    for (stampColIdx = -hwKSStamp; stampColIdx <= hwKSStamp; stampColIdx++) {
        for (stampRowIdx = -hwKSStamp; stampRowIdx <= hwKSStamp; stampRowIdx++) {
            pixelIdx  = stampColIdx + hwKSStamp + fwKSStamp * (stampRowIdx + hwKSStamp);
            q += vec[ncomp1][pixelIdx] * image[stampColIdx+xi+rPixX*(stampRowIdx+yi)];
        }
    }
    stamp->scprod[ncomp1+1] = q;

    return;
}

/**
 * @brief Perform an independent per-stamp kernel fit, sigma-clip outlier stamps
 *        by kernel sum, optionally do a global fit, and return a figure of merit.
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
double check_stamps(stamp_struct *stamps, int nS, float *imRef, float *imNoise) {

    int    nComps,stampIdx,matrixRowIdx,matrixColIdx,meritCount1,meritCount2,meritCount3;
    double sum=0,kmean,kstdev;
    double merit1,merit2,merit3,sig1,sig2,sig3;
    float *meritValues1,*meritValues2,*meritValues3,*ks;
    int    substampCenterX, substampCenterY, nks;

    double **matrix;
    int mat_size;
    int ncomp1, ncomp2, ncomp, nbg_vec;

    int ntestStamps;
    double       *testKerSol = NULL;
    stamp_struct *testStamps = NULL;

    /* kernel sum */
    ks  = (float *)calloc(nS, sizeof(float));
    nks = 0;

    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    mat_size   = ncomp1 * ncomp2 + nbg_vec + 1;

    if (verbose>=2) fprintf(stderr, " Mat_size0: %i ncomp2: %i ncomp1: %i nbg_vec: %i \n"
                            , mat_size, ncomp2, ncomp1, nbg_vec);

    /* for inital fit */
    nComps      = nCompKer + 1;

    for (stampIdx = 0; stampIdx < nS; stampIdx++) {

        substampCenterX    = stamps[stampIdx].xss[stamps[stampIdx].sscnt];
        substampCenterY    = stamps[stampIdx].yss[stamps[stampIdx].sscnt];

        /* extract check_mat to solve one particular stamp */
        for (matrixRowIdx = 1; matrixRowIdx <= nComps; matrixRowIdx++) {
            check_vec[matrixRowIdx] = stamps[stampIdx].scprod[matrixRowIdx];

            for (matrixColIdx = 1; matrixColIdx <= matrixRowIdx; matrixColIdx++) {
                check_mat[matrixRowIdx][matrixColIdx] = stamps[stampIdx].mat[matrixRowIdx][matrixColIdx];
                check_mat[matrixColIdx][matrixRowIdx] = check_mat[matrixRowIdx][matrixColIdx];
            }
        }

        /* fit stamp, the constant kernel coefficients end up in check_vec */
        solve_spd(check_mat, nComps, check_vec);

        /* find kernel sum */
        sum = check_vec[1];
        check_stack[stampIdx] = sum;
        stamps[stampIdx].norm = sum;
        ks[nks++]      = sum;

        if (verbose >= 2) fprintf(stderr, "    # %d    xss: %4i yss: %4i  ksum: %f\n", stampIdx,
                                  stamps[stampIdx].xss[stamps[stampIdx].sscnt],
                                  stamps[stampIdx].yss[stamps[stampIdx].sscnt], sum);
    }

    sigma_clip(ks, nks, &kmean, &kstdev, 10);

    fprintf(stderr, "    %.1f sigma clipped mean ksum : %.3f, stdev : %.3f, n : %i\n",
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
    if ((strncmp(forceConvolve, "b", 1)==0)) {

        /* allocate fitting matrix */
        matrix = (double **)calloc((mat_size + 1), sizeof(double *));
        for (stampIdx = 0; stampIdx <= mat_size; stampIdx++)
            matrix[stampIdx] = (double *)calloc((mat_size + 1), sizeof(double));

        /* allocate weight matrix */
        wxy = (double **)calloc(nS, sizeof(double *));
        for (stampIdx = 0; stampIdx < nS; stampIdx++)
            wxy[stampIdx] = (double *)calloc(ncomp2, sizeof(double));

        /* first find out how many good stamps to allocate */
        ntestStamps = 0;
        for (stampIdx = 0; stampIdx < nS; stampIdx++)
            if (stamps[stampIdx].diff < kerSigReject) {
                ntestStamps++;
            }
            else {
                if (verbose >= 2) fprintf(stderr, "    # %d    skipping xss: %4i yss: %4i ksum: %f sigma: %f\n", stampIdx,
                                          stamps[stampIdx].xss[stamps[stampIdx].sscnt],
                                          stamps[stampIdx].yss[stamps[stampIdx].sscnt],
                                          stamps[stampIdx].norm, stamps[stampIdx].diff);
            }

        /* then allocate test stamp structure */
        if(!(testStamps = (stamp_struct *)calloc(ntestStamps, sizeof(stamp_struct)))) {
            printf("Cannot Allocate Test Stamp List\n");
            exit (1);
        }
        testKerSol = (double *)calloc((nCompTotal+1), sizeof(double));

        /* and point test stamp structure to good stamps */
        ntestStamps = 0;
        for (stampIdx = 0; stampIdx < nS; stampIdx++)
            if (stamps[stampIdx].diff < kerSigReject)
                testStamps[ntestStamps++] = stamps[stampIdx];

        /* finally do fit */
        if (verbose >= 2) fprintf(stderr, " Expanding Test Matrix For Fit\n");
        build_matrix(testStamps, ntestStamps, matrix);
        build_scprod(testStamps, ntestStamps, imRef, testKerSol);
        solve_spd(matrix, mat_size, testKerSol);

        /* get the kernel sum to normalize figures of merit! */
        kmean = make_kernel(0, 0, testKerSol);

        /* determine figure of merit from good stamps */

        /* average of sum (diff**2 / value), ~variance */
        meritValues1 = (float *)calloc(ntestStamps, sizeof(float));

        /* standard deviation of pixel distribution */
        meritValues2 = (float *)calloc(ntestStamps, sizeof(float));

        /* noise sd based on histogram distribution width */
        meritValues3 = (float *)calloc(ntestStamps, sizeof(float));

        meritCount1 = 0;
        meritCount2 = 0;
        meritCount3 = 0;
        for (stampIdx = 0; stampIdx < ntestStamps; stampIdx++) {

            getStampSig(&testStamps[stampIdx], testKerSol, imNoise, &sig1, &sig2, &sig3);

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
        if (testKerSol)	 free(testKerSol);
        if (testStamps)	 free(testStamps);
        for (stampIdx = 0; stampIdx <= mat_size; stampIdx++)
            free(matrix[stampIdx]);
        for (stampIdx = 0; stampIdx < nS; stampIdx++)
            free(wxy[stampIdx]);
        free(matrix);
        free(wxy);

        free(meritValues1);
        free(meritValues2);
        free(meritValues3);
        free(ks);

        /* average value of figures of merit across stamps */
        fprintf(stderr, "    <var_merit> = %.3f, <sd_merit> = %.3f, <hist_merit> = %.3f\n", merit1, merit2, merit3);

        /* return what is asked for if possible, if not use backup */
        if (strncmp(figMerit, "v", 1)==0) {
            if (meritCount1 > 0) {
                return merit1;
            }
            else if (meritCount2 > 0) {
                return merit2;
            }
            else if (meritCount3 > 0) {
                return merit3;
            }
            else {
                return 666;
            }
        }
        else if (strncmp(figMerit, "s", 1)==0) {
            if (meritCount2 > 0) {
                return merit2;
            }
            else if (meritCount1 > 0) {
                return merit1;
            }
            else if (meritCount3 > 0) {
                return merit3;
            }
            else {
                return 666;
            }
        }
        else if (strncmp(figMerit, "h", 1)==0) {
            if (meritCount3 > 0) {
                return merit3;
            }
            else if (meritCount1 > 0) {
                return merit1;
            }
            else if (meritCount2 > 0) {
                return merit2;
            }
            else {
                return 666;
            }
        }
    }
    else
        return 0;

    return 0;
}

/**
 * @brief Alternative implementation of build_matrix() using a slightly
 *        different spatial coordinate normalisation.
 *
 * @details Functionally equivalent to build_matrix() but computes the
 * spatial polynomial weights wxy[istamp][k] using coordinates normalised as
 * fx = (xstamp - rPixX/2) / (rPixX/2) rather than the half-pixel-shifted
 * convention used in build_matrix().  This version is not called in the
 * default code path; build_matrix() is used instead.  Kept for reference and
 * potential future use.
 *
 * See build_matrix() for a full description of the matrix structure.
 *
 * @param stamps  Array of nS stamps with pre-built per-stamp matrices.
 * @param nS      Number of stamps.
 * @param matrix  Pre-allocated (mat_size+1)×(mat_size+1) output matrix,
 *                zeroed on entry; filled with the assembled normal-equations
 *                matrix on return.
 */
void build_matrix_new(stamp_struct *stamps, int nS, double **matrix) {
    
    int       mat_size,i,j,pixStamp,istamp,k,i1,i2,j1,j2,ibg,jbg,ivecbg,jj;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    **matrix0,p0,q,fx,fy;
    double    **vec;
    
    int ideg1, ideg2, xstamp, ystamp;
    double a1, a2;
    
    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    pixStamp   = fwKSStamp * fwKSStamp;
    
    mat_size = ncomp1 * ncomp2 + nbg_vec + 1;
    if (verbose >= 2) fprintf(stderr, " Mat_size: %i ncomp2: %i ncomp1: %i nbg_vec: %i \n",mat_size,ncomp2,ncomp1,nbg_vec);
    
    for(i = 0; i <= mat_size; i++)
        for(j = 0; j <= mat_size; j++)
            matrix[i][j] = 0.0;
    
    for(i = 0; i < nS; i++)
        for(j = 0; j < ncomp2; j++)
            wxy[i][j] = 0.0;
    
    for (istamp = 0; istamp < nS; istamp++) {
        /* skip over any bad stamps along the way */
        while(stamps[istamp].sscnt >= stamps[istamp].nss) { 
            ++istamp; 
            if (istamp >= nS) break;
        }
        if (istamp >= nS) break;
        
        vec    = stamps[istamp].vectors;
        xstamp = stamps[istamp].xss[stamps[istamp].sscnt];
        ystamp = stamps[istamp].yss[stamps[istamp].sscnt];
        
        /* build weight function *HERE*, implicitly giving bad stamps zero weight */
        /*   because we skip them over above... */
        k  = 0;
        
        fx = (xstamp - rPixX/2) / rPixX/2; 
        fy = (ystamp - rPixY/2) / rPixY/2; 
        
        for (ideg1 = 0, a1 = 1.0; ideg1 <= kerOrder; ideg1++, a1 *= fx)
            for (ideg2 = 0, a2 = 1.0; ideg2 <= kerOrder - ideg1; ideg2++, a2 *= fy)
                wxy[istamp][k++] = a1 * a2;

        matrix0 = stamps[istamp].mat;
        for (i = 0; i < ncomp; i++) {
            /* Decompose linear index i into (i1, i2):
               i1 = which Gaussian component (0 to ncomp1-1)
               i2 = which polynomial term within that component (0 to ncomp2-1) */
            i1 = i / ncomp2;
            i2 = i - i1 * ncomp2;

            for (j = 0; j <= i; j++) {
                /* Same decomposition for column index j */
                j1 = j / ncomp2;
                j2 = j - j1 * ncomp2;

                /* spatially weighted W_m1 and W_m2 integrals */
                matrix[i+2][j+2] += wxy[istamp][i2] * wxy[istamp][j2] * matrix0[i1+2][j1+2];
            }
        }
        
        matrix[1][1] += matrix0[1][1];
        for (i = 0; i < ncomp; i++) {
            i1 = i / ncomp2;
            i2 = i - i1 * ncomp2;
            matrix[i+2][1] += wxy[istamp][i2] * matrix0[i1+2][1];
        }
        
        for (ibg = 0; ibg < nbg_vec; ibg++) {
            i = ncomp + ibg + 1;
            ivecbg = ncomp1 + ibg + 1;
            for (i1 = 1; i1 < ncomp1 + 1; i1++) { 
                p0 = 0.0;
                
                /* integrate convolved images over all order backgrounds */
                for (k = 0; k < pixStamp; k++)
                    p0 += vec[i1][k] * vec[ivecbg][k];
                
                /* spatially weighted image * background terms */
                for (i2 = 0; i2 < ncomp2; i2++) {
                    jj = (i1 - 1) * ncomp2 + i2 + 1;
                    matrix[i+1][jj+1] += p0 * wxy[istamp][i2];
                }
            }
            
            p0 = 0.0;
            for (k = 0; k < pixStamp; k++)
                p0 += vec[0][k] * vec[ivecbg][k];
            matrix[i+1][1] += p0;
            
            /* background * background */
            for (jbg = 0;jbg <= ibg; jbg++) {
                for (k = 0, q = 0.0; k < pixStamp; k++)
                    q += vec[ivecbg][k] * vec[ncomp1+jbg+1][k];
                matrix[i+1][ncomp+jbg+2] += q;
            }
        }
    }
    
    /* fill lower half of matrix */
    for (i = 0; i < mat_size; i++) {
        for (j = 0; j <= i; j++) {
            matrix[j+1][i+1] = matrix[i+1][j+1];
            /* fprintf(stderr, "matrix[%i][%i]: %lf\n", i,j,matrix[i+1][j+1]); */
        }
    }
    
    return;
}

/**
 * @brief Assemble the global normal-equations matrix M for the spatially-varying
 *        kernel fit by accumulating spatially-weighted per-stamp cross-product
 *        matrices.
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
void build_matrix(stamp_struct *stamps, int nS, double **matrix) {

    int       mat_size,matrixRowIdx,matrixColIdx,pixStamp,stampIdx,polyTermIdx,gaussianCompIdx1,polyTermWithinGaussianIdx1,gaussianCompIdx2,polyTermIdx2,polyTermWithinGaussianIdx2,bgTermIdx,bgTermColIdx,bgVectorIdx,matrixColCalcIdx;
    int       polyTermIdx1;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    **matrix0,p0,q;
    double    **vec;
    float     rPixX2, rPixY2;

    int polyDegX, polyDegY, stampCenterX, stampCenterY;
    double polyBasisX, polyBasisY, normalizedX, normalizedY;

    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;

    pixStamp = fwKSStamp * fwKSStamp;
    rPixX2   = 0.5 * rPixX;
    rPixY2   = 0.5 * rPixY;

    mat_size = ncomp1 * ncomp2 + nbg_vec + 1;
    if (verbose >= 2) fprintf(stderr, " Mat_size: %i ncomp2: %i ncomp1: %i nbg_vec: %i \n",mat_size,ncomp2,ncomp1,nbg_vec);

    for(matrixRowIdx = 0; matrixRowIdx <= mat_size; matrixRowIdx++)
        for(matrixColIdx = 0; matrixColIdx <= mat_size; matrixColIdx++)
            matrix[matrixRowIdx][matrixColIdx] = 0.0;

    for(stampIdx = 0; stampIdx < nS; stampIdx++)
        for(polyTermIdx = 0; polyTermIdx < ncomp2; polyTermIdx++)
            wxy[stampIdx][polyTermIdx] = 0.0;

    for (stampIdx = 0; stampIdx < nS; stampIdx++) {
        /* skip over any bad stamps along the way */
        while(stamps[stampIdx].sscnt >= stamps[stampIdx].nss) {
            ++stampIdx;
            if (stampIdx >= nS) break;
        }
        if (stampIdx >= nS) break;

        vec    = stamps[stampIdx].vectors;
        stampCenterX = stamps[stampIdx].xss[stamps[stampIdx].sscnt];
        stampCenterY = stamps[stampIdx].yss[stamps[stampIdx].sscnt];
        /* RANGE FROM -1 to 1 */
        normalizedX     = (stampCenterX - rPixX2) / rPixX2;
        normalizedY     = (stampCenterY - rPixY2) / rPixY2;

        /* build weight function *HERE* */
        polyTermIdx  = 0;
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
            /* Decompose linear index polyTermIdx1 into (gaussianCompIdx1, polyTermWithinGaussianIdx1):
               gaussianCompIdx1 = which Gaussian component (0 to ncomp1-1)
               polyTermWithinGaussianIdx1 = which polynomial term within that component (0 to ncomp2-1) */
            gaussianCompIdx1 = polyTermIdx1 / ncomp2;
            polyTermWithinGaussianIdx1 = polyTermIdx1 - gaussianCompIdx1 * ncomp2;

            for (polyTermIdx2 = 0; polyTermIdx2 <= polyTermIdx1; polyTermIdx2++) {
                /* Same decomposition for column index polyTermIdx2 */
                gaussianCompIdx2 = polyTermIdx2 / ncomp2;
                polyTermWithinGaussianIdx2 = polyTermIdx2 - gaussianCompIdx2 * ncomp2;

                /* spatially weighted W_m1 and W_m2 integrals */
                matrix[polyTermWithinGaussianIdx1+2][polyTermWithinGaussianIdx2+2] += wxy[stampIdx][polyTermWithinGaussianIdx1] * wxy[stampIdx][polyTermWithinGaussianIdx2] * matrix0[gaussianCompIdx1+2][gaussianCompIdx2+2];
            }
        }

        matrix[1][1] += matrix0[1][1];
        for (polyTermIdx1 = 0; polyTermIdx1 < ncomp; polyTermIdx1++) {
            /* Index decomposition for background terms */
            gaussianCompIdx1 = polyTermIdx1 / ncomp2;
            polyTermWithinGaussianIdx1 = polyTermIdx1 - gaussianCompIdx1 * ncomp2;
            matrix[polyTermWithinGaussianIdx1+2][1] += wxy[stampIdx][polyTermWithinGaussianIdx1] * matrix0[gaussianCompIdx1+2][1];
        }

        for (bgTermIdx = 0; bgTermIdx < nbg_vec; bgTermIdx++) {
            polyTermIdx1 = ncomp + bgTermIdx + 1;
            bgVectorIdx = ncomp1 + bgTermIdx + 1;
            for (gaussianCompIdx1 = 1; gaussianCompIdx1 < ncomp1 + 1; gaussianCompIdx1++) {
                p0 = 0.0;

                /* integrate convolved images over all order backgrounds */
                for (polyTermIdx = 0; polyTermIdx < pixStamp; polyTermIdx++)
                    p0 += vec[gaussianCompIdx1][polyTermIdx] * vec[bgVectorIdx][polyTermIdx];

                /* spatially weighted image * background terms */
                for (polyTermIdx2 = 0; polyTermIdx2 < ncomp2; polyTermIdx2++) {
                    matrixColCalcIdx = (gaussianCompIdx1 - 1) * ncomp2 + polyTermIdx2 + 1;
                    matrix[polyTermIdx1+1][matrixColCalcIdx+1] += p0 * wxy[stampIdx][polyTermIdx2];
                }
            }

            p0 = 0.0;
            for (polyTermIdx = 0; polyTermIdx < pixStamp; polyTermIdx++)
                p0 += vec[0][polyTermIdx] * vec[bgVectorIdx][polyTermIdx];
            matrix[polyTermIdx1+1][1] += p0;

            /* background * background */
            for (bgTermColIdx = 0;bgTermColIdx <= bgTermIdx; bgTermColIdx++) {
                for (polyTermIdx = 0, q = 0.0; polyTermIdx < pixStamp; polyTermIdx++)
                    q += vec[bgVectorIdx][polyTermIdx] * vec[ncomp1+bgTermColIdx+1][polyTermIdx];
                matrix[polyTermIdx1+1][ncomp+bgTermColIdx+2] += q;
            }
        }
    }

    /* fill lower half of matrix */
    for (matrixRowIdx = 0; matrixRowIdx < mat_size; matrixRowIdx++) {
        for (matrixColIdx = 0; matrixColIdx <= matrixRowIdx; matrixColIdx++) {
            matrix[matrixColIdx+1][matrixRowIdx+1] = matrix[matrixRowIdx+1][matrixColIdx+1];
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
 * build_matrix() fills M and this function fills b, ludcmp/lubksb solve M·a=b
 * to yield the kernel coefficient vector in @p kernelSol.
 *
 * @param stamps     Array of nS stamps.
 * @param nS         Number of stamps.
 * @param image      Full-frame reference image (used for the background terms).
 * @param kernelSol  Output: right-hand-side vector of length nCompTotal+1;
 *                   zeroed on entry, filled on return.
 */
void build_scprod(stamp_struct *stamps, int nS, float *image, double *kernelSol) {

    int       stampIdx,stampColIdx,stampRowIdx,stampCenterX,stampCenterY,gaussianCompIdx,polyTermIdx,pixelIdx,bgVecIdx,solIdx,recomputedColIdx;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    scalarProduct,innerProduct;
    double    **vectorArray;

    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;


    for (solIdx = 0; solIdx <= ncomp + nbg_vec + 1; solIdx++)
        kernelSol[solIdx]=0.0;

    for (stampIdx = 0; stampIdx < nS; stampIdx++) {
        /* skip over any bad stamps along the way */
        while(stamps[stampIdx].sscnt >= stamps[stampIdx].nss) {
            ++stampIdx;
            if (stampIdx >= nS) break;
        }
        if (stampIdx >= nS) break;

        vectorArray= stamps[stampIdx].vectors;
        stampCenterX = stamps[stampIdx].xss[stamps[stampIdx].sscnt];
        stampCenterY = stamps[stampIdx].yss[stamps[stampIdx].sscnt];

        scalarProduct = stamps[stampIdx].scprod[1];
        kernelSol[1] += scalarProduct;

        /* spatially weighted convolved image * ref image */
        for (gaussianCompIdx = 1; gaussianCompIdx < ncomp1 + 1; gaussianCompIdx++) {
            scalarProduct = stamps[stampIdx].scprod[gaussianCompIdx+1];
            for (polyTermIdx = 0; polyTermIdx < ncomp2; polyTermIdx++) {
                recomputedColIdx = (gaussianCompIdx-1) * ncomp2 + polyTermIdx + 1;
                /* no need for weighting here */
                kernelSol[recomputedColIdx+1] += scalarProduct * wxy[stampIdx][polyTermIdx];
            }
        }

        /* spatially weighted bg model convolved with ref image */
        for (bgVecIdx = 0; bgVecIdx < nbg_vec; bgVecIdx++) {
            innerProduct = 0.0;
            for (stampColIdx = -hwKSStamp; stampColIdx <= hwKSStamp; stampColIdx++) {
                for (stampRowIdx = -hwKSStamp; stampRowIdx <= hwKSStamp;stampRowIdx++) {
                    pixelIdx  = stampColIdx + hwKSStamp + fwKSStamp * (stampRowIdx + hwKSStamp);
                    innerProduct += vectorArray[ncomp1+bgVecIdx+1][pixelIdx] * image[stampColIdx+stampCenterX+rPixX*(stampRowIdx+stampCenterY)];

                }
            }
            kernelSol[ncomp+bgVecIdx+2] += innerProduct;
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
void getFinalStampSig(stamp_struct *stamp, float *imDiff, float *imNoise, double *sig) {
    int colIdx, rowIdx, pixelIdx, validPixelCount=0;
    int xPixel, xRegion = stamp->xss[stamp->sscnt];
    int yPixel, yRegion = stamp->yss[stamp->sscnt];
    float diffData, noiseInv;

    *sig = 0;

    for (rowIdx = 0; rowIdx < fwKSStamp; rowIdx++) {
        yPixel = yRegion - hwKSStamp + rowIdx;

        for (colIdx = 0; colIdx < fwKSStamp; colIdx++) {
            xPixel = xRegion - hwKSStamp + colIdx;

            pixelIdx   = xPixel+rPixX*yPixel;
            diffData  = imDiff[pixelIdx];
            noiseInv = 1. / imNoise[pixelIdx];

            /* this shouldn't be the case, but just in case... */
            if (mRData[pixelIdx] & FLAG_INPUT_ISBAD)
                continue;

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
void getStampSig(stamp_struct *stamp, double *kernelSol, float *imNoise, double *sig1, double *sig2, double *sig3) {
    int colIdx, rowIdx, stampPixelIdx, substampIdx, validPixelCount, xRegion, yRegion, xPixel, yPixel;
    double cSum, cMean, cMedian, cMode, cLfwhm;
    double *referenceImage, templateData, imageData, noiseData, diffValue, bgValue;

    /* info */
    substampIdx   = stamp->sscnt;
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

            stampPixelIdx  = colIdx+rowIdx*fwKSStamp;

            templateData = temp[stampPixelIdx];
            imageData = referenceImage[stampPixelIdx];
            noiseData = imNoise[xPixel+rPixX*yPixel];

            diffValue = templateData - imageData + bgValue;

            if ((mRData[xPixel+rPixX*yPixel] & FLAG_INPUT_ISBAD) || (fabs(imageData) <= ZEROVAL)) {
                continue;
            }
            else {
                temp[stampPixelIdx] = diffValue;
            }

            /* check for NaN */
            if ((templateData*0.0 != 0.0) || (imageData*0.0 != 0.0)) {
                mRData[xPixel+rPixX*yPixel] |= (FLAG_INPUT_ISBAD | FLAG_ISNAN);
                continue;
            }

            validPixelCount++;
            *sig1 += diffValue * diffValue / noiseData;
            /*fprintf(stderr, "OK %d %d : %f %f %f\n", xPixel, yPixel, templateData, imageData, noiseData);*/
        }
    }
    if (validPixelCount > 0) {
        *sig1 /= validPixelCount;
        if (*sig1 >= MAXVAL)
            *sig1 = -1;
    }
    else
        *sig1 = -1;


    /* don't do think unless you need to! */
    if (strncmp(figMerit, "v", 1)!=0) {
        if (getStampStats3(temp, xRegion - hwKSStamp, yRegion - hwKSStamp, fwKSStamp, fwKSStamp,
                           &cSum, &cMean, &cMedian, &cMode, sig2, sig3, &cLfwhm, 0x0, 0xffff, 5)) {
            *sig2 = -1;
            *sig3 = -1;
        }
        else if (*sig2 < 0 || *sig2 >= MAXVAL)
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
char check_again(stamp_struct *stamps, double *kernelSol, float *imConv, float *imRef, float *imNoise, double *meansigSubstamps, double *scatterSubstamps, int *NskippedSubstamps) {
    
    int    istamp,nss,scnt;
    double sig,mean,stdev;
    char   check;
    double sig1, sig2, sig3;
    float  *ss;
    
    ss  = (float *)calloc(nS, sizeof(float));
    nss = 0;
    
    sig = 0;
    check = 0;
    mean = stdev = 0.0;
    *NskippedSubstamps=0;
    
    for (istamp = 0; istamp < nS; istamp++) {
        
        /* if was fit with a good legit substamp */
        if (stamps[istamp].sscnt < stamps[istamp].nss) {
            
            getStampSig(&stamps[istamp], kernelSol, imNoise, &sig1, &sig2, &sig3);
            
            if ((strncmp(figMerit, "v", 1)==0 && (sig1 == -1)) ||
                (strncmp(figMerit, "s", 1)==0 && (sig2 == -1)) ||
                (strncmp(figMerit, "h", 1)==0 && (sig3 == -1))) {
                
                /* something went wrong with this one... */
                if (verbose>=2) fprintf(stderr, "\n    # %d    xss: %4i yss: %4i sig: %6.3f sscnt: %2i nss: %2i ITERATE substamp (BAD)\n",
                                        istamp,
                                        stamps[istamp].xss[stamps[istamp].sscnt],
                                        stamps[istamp].yss[stamps[istamp].sscnt],
                                        sig, stamps[istamp].sscnt, stamps[istamp].nss);
                
                stamps[istamp].sscnt++;
                fillStamp(&stamps[istamp], imConv, imRef);
                if (verbose>=2) fprintf(stderr, "\n");
                
                check = 1;
                
            } else {
                if (strncmp(figMerit, "v", 1)==0)
                    sig = sig1;
                else if (strncmp(figMerit, "s", 1)==0)
                    sig = sig2;
                else if (strncmp(figMerit, "h", 1)==0)
                    sig = sig3;
                
                if (verbose>=2) fprintf(stderr, "    # %d    xss: %4i yss: %4i sig: %6.3f sscnt: %2i nss: %2i OK\n",
                                        istamp,
                                        stamps[istamp].xss[stamps[istamp].sscnt],
                                        stamps[istamp].yss[stamps[istamp].sscnt],
                                        sig, stamps[istamp].sscnt, stamps[istamp].nss);
                
                stamps[istamp].chi2 = sig;
                ss[nss++]           = sig;
                
            }
        } else {
            (*NskippedSubstamps)++;
            if (verbose>=2) fprintf(stderr, "    xs : %4i ys : %4i skipping... \n",stamps[istamp].x, stamps[istamp].y);
        }
    }
    
    sigma_clip(ss, nss, &mean, &stdev, 10);
    fprintf(stderr, "    Mean sig: %6.3f stdev: %6.3f\n", mean, stdev);
    fprintf(stderr, "    Iterating through stamps with sig > %.3f\n", mean + kerSigReject * stdev);
    
    /* save the mean and scatter so that it can be saved in the fits header */
    (*meansigSubstamps)=mean;
    (*scatterSubstamps)=stdev;
    
    scnt = 0;
    for (istamp = 0; istamp < nS; istamp++) {
        /* if currently represented by a good substamp */
        if (stamps[istamp].sscnt < stamps[istamp].nss) {
            
            /* no fabs() here, keep good stamps kerSigReject on the low side! */
            if ((stamps[istamp].chi2 - mean) > kerSigReject * stdev) { 
                if (verbose>=2) fprintf(stderr, "\n    # %d    xss: %4i yss: %4i sig: %6.3f sscnt: %2i nss: %2i ITERATE substamp (poor sig)\n",
                                        istamp,
                                        stamps[istamp].xss[stamps[istamp].sscnt],
                                        stamps[istamp].yss[stamps[istamp].sscnt],
                                        stamps[istamp].chi2,
                                        stamps[istamp].sscnt,stamps[istamp].nss);
                
                stamps[istamp].sscnt++;
                scnt += (!(fillStamp(&stamps[istamp], imConv, imRef)));
                if (verbose>=2) fprintf(stderr, "\n");
                
                check = 1;
            }
            else
                scnt += 1;
        }
    }
    
    fprintf(stderr, "    %d out of %d stamps remain\n", scnt, nS);
    
    free(ss);
    return check;
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
static double make_kernel_local(int xi, int yi, double *kernelSol,
                                double *lkernel, double *lkernel_coeffs) {
    int    i1, k, ix, iy, i;
    double ax, ay, sum_kernel;
    double xf, yf;

    k  = 2;
    xf = (xi - 0.5 * rPixX) / (0.5 * rPixX);
    yf = (yi - 0.5 * rPixY) / (0.5 * rPixY);

    for (i1 = 1; i1 < nCompKer; i1++) {
        lkernel_coeffs[i1] = 0.0;
        ax = 1.0;
        for (ix = 0; ix <= kerOrder; ix++) {
            ay = 1.0;
            for (iy = 0; iy <= kerOrder - ix; iy++) {
                lkernel_coeffs[i1] += kernelSol[k++] * ax * ay;
                ay *= yf;
            }
            ax *= xf;
        }
    }
    lkernel_coeffs[0] = kernelSol[1];

    for (i = 0; i < fwKernel * fwKernel; i++)
        lkernel[i] = 0.0;

    sum_kernel = 0.0;
    for (i = 0; i < fwKernel * fwKernel; i++) {
        for (i1 = 0; i1 < nCompKer; i1++)
            lkernel[i] += lkernel_coeffs[i1] * kernel_vec[i1][i];
        sum_kernel += lkernel[i];
    }
    return sum_kernel;
}

#ifdef USE_FFTW
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
 *
 * FFTW plan creation is serialised via #pragma omp critical(fftw_plan);
 * plan execution and the pixel-sum loop are parallelised with OpenMP.
 * CFITSIO is not called here; all FITS I/O occurs in main.c.
 *
 * Pixel values are produced by the FFT path.  Mask propagation and variance
 * are computed by the same block-centre direct loop used in spatial_convolve(),
 * with the pixel accumulation (q +=) removed since cRdata is already filled.
 *
 * @param image      Full-frame input image (region pixel buffer).
 * @param variance   Pointer to variance image; updated in-place if non-NULL.
 * @param xSize      Region width in pixels.
 * @param ySize      Region height in pixels.
 * @param kernelSol  Kernel solution vector from fitKernel().
 * @param cRdata     Output convolved image (pre-allocated by caller).
 * @param cMask      Input pixel mask for mask propagation.
 */
static void spatial_convolve_fft(float *image, float **variance,
                                  int xSize, int ySize,
                                  double *kernelSol,
                                  float *cRdata, int *cMask) {
    /* --- variable declarations (all at top for safe goto use) --- */
    int    ncomp2, nEffConv, fft_nx, fft_ny, nc_fft;
    int    i, j, p, ik, jk, ii, jj, xp, yp, ni, ix, iy;
    int    i1, i2, j2, i0, j0, j1, nsteps_x_loc, nsteps_y_loc;
    int    ic, jc, nc, mbit, dovar, cidx;
    double ax, ay, out, xf, yf, kk, aks, uks, qv, norm, rPixX2, rPixY2;
    double *effKernel, *real_buf, *lkernel, *lkernel_coeffs;
    fftw_complex *img_fft, *ker_fft;
    fftw_plan plan_fwd, plan_inv;
    float **effConv, *vData;

    effKernel     = NULL; real_buf = NULL;
    lkernel       = NULL; lkernel_coeffs = NULL;
    img_fft       = NULL; ker_fft  = NULL;
    plan_fwd      = NULL; plan_inv = NULL;
    effConv       = NULL; vData    = NULL;

    ncomp2   = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    nEffConv = 1 + ncomp2;
    dovar    = (*variance != NULL);

    if (dovar) {
        vData = (float *)calloc((size_t)xSize * ySize, sizeof(float));
        if (!vData) {
            fprintf(stderr, "spatial_convolve_fft: out of memory (vData)\n");
            return;
        }
    }

    /* --- FFT dimensions: pad to avoid circular wrap-around --- */
    fft_nx = next_fftw_size(xSize + fwKernel - 1);
    fft_ny = next_fftw_size(ySize + fwKernel - 1);
    nc_fft = fft_nx * (fft_ny / 2 + 1);

    real_buf  = (double *)fftw_malloc((size_t)fft_nx * fft_ny * sizeof(double));
    img_fft   = (fftw_complex *)fftw_malloc((size_t)nc_fft * sizeof(fftw_complex));
    ker_fft   = (fftw_complex *)fftw_malloc((size_t)nc_fft * sizeof(fftw_complex));
    effKernel = (double *)malloc((size_t)fwKernel * fwKernel * sizeof(double));

    if (!real_buf || !img_fft || !ker_fft || !effKernel) {
        fprintf(stderr, "spatial_convolve_fft: out of memory (FFT buffers)\n");
        goto cleanup_fft;
    }

    effConv = (float **)calloc((size_t)nEffConv, sizeof(float *));
    if (!effConv) {
        fprintf(stderr, "spatial_convolve_fft: out of memory (effConv)\n");
        goto cleanup_fft;
    }
    for (i = 0; i < nEffConv; i++) {
        effConv[i] = (float *)malloc((size_t)xSize * ySize * sizeof(float));
        if (!effConv[i]) {
            fprintf(stderr, "spatial_convolve_fft: out of memory (effConv[%d])\n", i);
            goto cleanup_fft;
        }
    }

    /* --- Create FFTW plans (not thread-safe: serialise with critical) --- */
#ifdef _OPENMP
    #pragma omp critical(fftw_plan)
#endif
    {
        /* FFTW_ESTIMATE avoids timing runs, keeping plan creation fast. */
        plan_fwd = fftw_plan_dft_r2c_2d(fft_ny, fft_nx,
                                         real_buf, ker_fft, FFTW_ESTIMATE);
        plan_inv = fftw_plan_dft_c2r_2d(fft_ny, fft_nx,
                                         ker_fft, real_buf, FFTW_ESTIMATE);
    }
    if (!plan_fwd || !plan_inv) {
        fprintf(stderr, "spatial_convolve_fft: FFTW plan creation failed\n");
        goto cleanup_fft;
    }

    /* --- FFT the image once; save spectrum in img_fft --- */
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
            double *kd = (double *)ker_fft;
            const double *id = (const double *)img_fft;
            double reval, imval;
            for (p = 0; p < nc_fft; p++) {
                reval       = kd[2*p]*id[2*p]   - kd[2*p+1]*id[2*p+1];
                imval       = kd[2*p]*id[2*p+1] + kd[2*p+1]*id[2*p];
                kd[2*p]     = reval;
                kd[2*p+1]   = imval;
            }
        }
        fftw_execute_dft_c2r(plan_inv, ker_fft, real_buf);

        /* Normalise (FFTW does not normalise) and store valid region */
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
        fftw_destroy_plan(plan_fwd); plan_fwd = NULL;
        fftw_destroy_plan(plan_inv); plan_inv = NULL;
    }
    fftw_free(real_buf); real_buf = NULL;
    fftw_free(img_fft);  img_fft  = NULL;
    fftw_free(ker_fft);  ker_fft  = NULL;
    free(effKernel);     effKernel = NULL;

    /* --- Pixel-wise weighted sum (OpenMP parallel over output pixels) ---
     * output(x,y) = effConv[0][ni]
     *             + Σⱼ Pⱼ(xf,yf) * effConv[j+1][ni]
     * where Pⱼ = xf^ix * yf^iy in the same (ix,iy) order as make_kernel_local.
     * Each ni = xp + xSize*yp is unique; no data race on cRdata or mRData. */
    rPixX2 = 0.5 * xSize;
    rPixY2 = 0.5 * ySize;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) collapse(2) \
        private(ni, xf, yf, out, ax, ay, j, ix, iy)
#endif
    for (yp = hwKernel; yp < ySize - hwKernel; yp++) {
        for (xp = hwKernel; xp < xSize - hwKernel; xp++) {
            ni  = xp + xSize * yp;
            xf  = (xp - rPixX2) / rPixX2;
            yf  = (yp - rPixY2) / rPixY2;

            out = effConv[0][ni];
            j   = 0;
            ax  = 1.0;
            for (ix = 0; ix <= kerOrder; ix++) {
                ay = 1.0;
                for (iy = 0; iy <= kerOrder - ix; iy++) {
                    out += ax * ay * effConv[j + 1][ni];
                    ay  *= yf;
                    j++;
                }
                ax *= xf;
            }
            cRdata[ni] = (float)out;
        }
    }

    /* --- Free effective convolution images --- */
    for (i = 0; i < nEffConv; i++) { free(effConv[i]); effConv[i] = NULL; }
    free(effConv); effConv = NULL;

    /* --- Mask propagation + variance: block-centre direct loop ---
     * Structurally identical to the non-OMP path in spatial_convolve() with
     * the pixel accumulation (q +=) removed; cRdata is already filled above.
     * Variance is propagated exactly as in the original. */
    lkernel        = (double *)calloc((size_t)fwKernel * fwKernel, sizeof(double));
    lkernel_coeffs = (double *)calloc((size_t)nCompKer,             sizeof(double));
    if (!lkernel || !lkernel_coeffs) {
        fprintf(stderr, "spatial_convolve_fft: out of memory (mask loop)\n");
        goto cleanup_fft;
    }

    nsteps_x_loc = (int)ceil((double)xSize / (double)kcStep);
    nsteps_y_loc = (int)ceil((double)ySize / (double)kcStep);

    for (j1 = 0; j1 < nsteps_y_loc; j1++) {
        j0 = j1 * kcStep + hwKernel;

        for (i = 0; i < nsteps_x_loc; i++) {
            i0 = i * kcStep + hwKernel;

            make_kernel_local(i0 + hwKernel, j0 + hwKernel, kernelSol,
                              lkernel, lkernel_coeffs);

            for (j2 = 0; j2 < kcStep; j2++) {
                j = j0 + j2;
                if (j >= ySize - hwKernel) break;

                for (i2 = 0; i2 < kcStep; i2++) {
                    i1 = i0 + i2;
                    if (i1 >= xSize - hwKernel) break;

                    ni   = i1 + xSize * j;
                    qv   = aks = uks = 0.0;
                    mbit = 0x0;

                    for (jc = j - hwKernel; jc <= j + hwKernel; jc++) {
                        jk = j - jc + hwKernel;
                        for (ic = i1 - hwKernel; ic <= i1 + hwKernel; ic++) {
                            ik = i1 - ic + hwKernel;
                            nc = ic + xSize * jc;
                            kk = lkernel[ik + jk * fwKernel];

                            if (dovar) {
                                if (convolveVariance)
                                    qv += (*variance)[nc] * kk;
                                else
                                    qv += (*variance)[nc] * kk * kk;
                            }

                            mbit |= cMask[nc];
                            aks  += fabs(kk);
                            if (!(cMask[nc] & FLAG_INPUT_ISBAD))
                                uks += fabs(kk);
                        }
                    }

                    if (dovar) vData[ni] = (float)qv;

                    mRData[ni] |= cMask[ni];
                    mRData[ni] |= FLAG_OUTPUT_ISBAD * ((cMask[ni] & FLAG_INPUT_ISBAD) > 0);

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

    free(lkernel);        lkernel        = NULL;
    free(lkernel_coeffs); lkernel_coeffs = NULL;

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
        if (plan_fwd) { fftw_destroy_plan(plan_fwd); plan_fwd = NULL; }
        if (plan_inv) { fftw_destroy_plan(plan_inv); plan_inv = NULL; }
    }
    if (real_buf)  { fftw_free(real_buf);  real_buf  = NULL; }
    if (img_fft)   { fftw_free(img_fft);   img_fft   = NULL; }
    if (ker_fft)   { fftw_free(ker_fft);   ker_fft   = NULL; }
    if (effKernel) { free(effKernel);      effKernel = NULL; }
    if (effConv) {
        for (i = 0; i < nEffConv; i++) free(effConv[i]);
        free(effConv); effConv = NULL;
    }
    if (lkernel)        { free(lkernel);        lkernel        = NULL; }
    if (lkernel_coeffs) { free(lkernel_coeffs); lkernel_coeffs = NULL; }
    free(vData); vData = NULL;
}
#endif /* USE_FFTW */

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
 * - FFT acceleration: O(N log N) precomputed basis convolutions (preferred for large images)
 *
 * Dispatches to spatial_convolve_fft() if USE_FFTW is enabled and FFTW3 is available;
 * otherwise uses direct O(N·k²) implementation.
 *
 * @param[in] image      Flat image array (e.g., template) to be convolved.
 * @param[in,out] variance Optional variance map (pointer-to-pointer); if non-NULL,
 *                        noise propagation is computed in-place.
 * @param[in] xSize      Image width in pixels.
 * @param[in] ySize      Image height in pixels.
 * @param[in] kernelSol  Kernel solution vector cᵢ from fitKernel(); dimensionality nCompKer.
 * @param[out] cRdata    Output convolved image (pre-allocated by caller; xSize×ySize).
 * @param[in] cMask      Input pixel mask for mask propagation (0 = bad, non-zero = good).
 *
 * @note Calls make_kernel() to evaluate K(x,y) at each pixel; expensive operation.
 *       For multi-threaded execution, loop indices and temporary buffers are thread-local.
 * @see spatial_convolve_fft() for FFT-accelerated variant
 * @see fitKernel() for kernel solution computation
 */
void spatial_convolve(float *image, float **variance, int xSize, int ySize, double *kernelSol, float *cRdata, int *cMask) {
#ifdef USE_FFTW
    spatial_convolve_fft(image, variance, xSize, ySize, kernelSol, cRdata, cMask);
    return;
#endif

    int       kernelStepColIdx,kernelStepRowIdx,pixelWithinStepColIdx,pixelWithinStepRowIdx,nsteps_x,nsteps_y,pixelX,pixelY,kernelStepOriginX,kernelStepOriginY,kernelCenterColIdx,kernelCenterRowIdx,kernelArrayColIdx,kernelArrayRowIdx,neighborPixelIdx,outputPixelIdx,maskedPixelFlags,hasVarianceImage;
    double    convolvedValue, convolvedVariance, kernelValue, absKernelSum, unmaskedKernelSum;
    float     *vData=NULL;

    if ((*variance) == NULL)
        hasVarianceImage = 0;
    else
        hasVarianceImage = 1;

    if (hasVarianceImage) {
        if ( !(vData = (float *)calloc(xSize*ySize, sizeof(float)))) {
            return;
        }
    }

    nsteps_x = ceil((double)(xSize)/(double)kcStep);
    nsteps_y = ceil((double)(ySize)/(double)kcStep);

#ifdef _OPENMP
    /* Each thread gets private copies of the loop indices and accumulators.
       kernelStepRowIdx is automatically private as the omp-for loop variable.
       lkernel / lkernel_coeffs are allocated per-thread inside the parallel
       block so that concurrent make_kernel_local() calls are race-free. */
    #pragma omp parallel \
        private(kernelStepColIdx, pixelWithinStepColIdx, pixelWithinStepRowIdx, pixelX, pixelY, kernelStepOriginX, kernelStepOriginY, kernelCenterColIdx, kernelCenterRowIdx, kernelArrayColIdx, kernelArrayRowIdx, neighborPixelIdx, outputPixelIdx, maskedPixelFlags, convolvedValue, convolvedVariance, kernelValue, absKernelSum, unmaskedKernelSum)
    {
        double *lkernel       = (double *)calloc(fwKernel * fwKernel, sizeof(double));
        double *lkernel_coeffs = (double *)calloc(nCompKer, sizeof(double));

        #pragma omp for schedule(dynamic)
        for (kernelStepRowIdx = 0; kernelStepRowIdx < nsteps_y; kernelStepRowIdx++) {
            kernelStepOriginY = kernelStepRowIdx * kcStep + hwKernel;

            for (kernelStepColIdx = 0; kernelStepColIdx < nsteps_x; kernelStepColIdx++) {
                kernelStepOriginX = kernelStepColIdx * kcStep + hwKernel;

                make_kernel_local(kernelStepOriginX + hwKernel, kernelStepOriginY + hwKernel, kernelSol,
                                  lkernel, lkernel_coeffs);

                for (pixelWithinStepRowIdx = 0; pixelWithinStepRowIdx < kcStep; pixelWithinStepRowIdx++) {
                    pixelY = kernelStepOriginY + pixelWithinStepRowIdx;
                    if (pixelY >= ySize - hwKernel) break;

                    for (pixelWithinStepColIdx = 0; pixelWithinStepColIdx < kcStep; pixelWithinStepColIdx++) {
                        pixelX = kernelStepOriginX + pixelWithinStepColIdx;
                        if (pixelX >= xSize - hwKernel) break;

                        outputPixelIdx = pixelX + xSize * pixelY;
                        convolvedVariance = convolvedValue = absKernelSum = unmaskedKernelSum = 0.0;
                        maskedPixelFlags = 0x0;
                        for (kernelCenterRowIdx = pixelY - hwKernel; kernelCenterRowIdx <= pixelY + hwKernel; kernelCenterRowIdx++) {
                            kernelArrayRowIdx = pixelY - kernelCenterRowIdx + hwKernel;

                            for (kernelCenterColIdx = pixelX - hwKernel; kernelCenterColIdx <= pixelX + hwKernel; kernelCenterColIdx++) {
                                kernelArrayColIdx = pixelX - kernelCenterColIdx + hwKernel;

                                neighborPixelIdx = kernelCenterColIdx + xSize * kernelCenterRowIdx;
                                kernelValue = lkernel[kernelArrayColIdx + kernelArrayRowIdx * fwKernel];

                                convolvedValue    += image[neighborPixelIdx] * kernelValue;
                                if (hasVarianceImage) {
                                    if (convolveVariance)
                                        convolvedVariance += (*variance)[neighborPixelIdx] * kernelValue;
                                    else
                                        convolvedVariance += (*variance)[neighborPixelIdx] * kernelValue * kernelValue;
                                }

                                maskedPixelFlags |= cMask[neighborPixelIdx];
                                absKernelSum  += fabs(kernelValue);
                                if (!(cMask[neighborPixelIdx] & FLAG_INPUT_ISBAD))
                                    unmaskedKernelSum += fabs(kernelValue);
                            }
                        }

                        cRdata[outputPixelIdx] = convolvedValue;
                        if (hasVarianceImage)
                            vData[outputPixelIdx] = convolvedVariance;

                        mRData[outputPixelIdx] |= cMask[outputPixelIdx];
                        mRData[outputPixelIdx] |= FLAG_OUTPUT_ISBAD * ((cMask[outputPixelIdx] & FLAG_INPUT_ISBAD) > 0);

                        if (maskedPixelFlags) {
                            if ((unmaskedKernelSum / absKernelSum) < kerFracMask)
                                mRData[outputPixelIdx] |= (FLAG_OUTPUT_ISBAD | FLAG_BAD_CONV);
                            else
                                mRData[outputPixelIdx] |= FLAG_OK_CONV;
                        }
                    }
                }
            }
        }

        free(lkernel);
        free(lkernel_coeffs);
    }
#else
    for (kernelStepRowIdx = 0; kernelStepRowIdx < nsteps_y; kernelStepRowIdx++) {
        kernelStepOriginY = kernelStepRowIdx * kcStep + hwKernel;

        for(kernelStepColIdx = 0; kernelStepColIdx < nsteps_x; kernelStepColIdx++) {
            kernelStepOriginX = kernelStepColIdx * kcStep + hwKernel;

            make_kernel(kernelStepOriginX + hwKernel, kernelStepOriginY + hwKernel, kernelSol);

            for (pixelWithinStepRowIdx = 0; pixelWithinStepRowIdx < kcStep; pixelWithinStepRowIdx++) {
                pixelY = kernelStepOriginY + pixelWithinStepRowIdx;
                if ( pixelY >= ySize - hwKernel) break;

                for (pixelWithinStepColIdx = 0; pixelWithinStepColIdx < kcStep; pixelWithinStepColIdx++) {
                    pixelX = kernelStepOriginX + pixelWithinStepColIdx;
                    if (pixelX >= xSize - hwKernel) break;

                    outputPixelIdx = pixelX+xSize*pixelY;
                    convolvedVariance = convolvedValue = absKernelSum = unmaskedKernelSum = 0.0;
                    maskedPixelFlags = 0x0;
                    for (kernelCenterRowIdx = pixelY - hwKernel; kernelCenterRowIdx <= pixelY + hwKernel; kernelCenterRowIdx++) {
                        kernelArrayRowIdx = pixelY - kernelCenterRowIdx + hwKernel;

                        for (kernelCenterColIdx = pixelX - hwKernel; kernelCenterColIdx <= pixelX + hwKernel; kernelCenterColIdx++) {
                            kernelArrayColIdx = pixelX - kernelCenterColIdx + hwKernel;

                            neighborPixelIdx = kernelCenterColIdx+xSize*kernelCenterRowIdx;
                            kernelValue = kernel[kernelArrayColIdx+kernelArrayRowIdx*fwKernel];

                            convolvedValue     += image[neighborPixelIdx] * kernelValue;
                            if (hasVarianceImage) {
                                if (convolveVariance)
                                    convolvedVariance += (*variance)[neighborPixelIdx] * kernelValue;
                                else
                                    convolvedVariance += (*variance)[neighborPixelIdx] * kernelValue * kernelValue;
                            }

                            maskedPixelFlags  |= cMask[neighborPixelIdx];
                            absKernelSum   += fabs(kernelValue);
                            if (!(cMask[neighborPixelIdx] & FLAG_INPUT_ISBAD)) {
                                unmaskedKernelSum += fabs(kernelValue);
                            }
                        }
                    }

                    cRdata[outputPixelIdx]   = convolvedValue;
                    if (hasVarianceImage)
                        vData[outputPixelIdx] = convolvedVariance;

                    /* mask propagation changed in 5.1.9 */
                    /* mRData[outputPixelIdx]  |= maskedPixelFlags; */
                    /* mRData[outputPixelIdx]  |= FLAG_OK_CONV      * (maskedPixelFlags > 0);*/
                    mRData[outputPixelIdx]  |= cMask[outputPixelIdx];
                    mRData[outputPixelIdx]  |= FLAG_OUTPUT_ISBAD * ((cMask[outputPixelIdx] & FLAG_INPUT_ISBAD) > 0);

                    if (maskedPixelFlags) {
                        if ((unmaskedKernelSum / absKernelSum) < kerFracMask) {
                            mRData[outputPixelIdx] |= (FLAG_OUTPUT_ISBAD | FLAG_BAD_CONV);
                        }
                        else {
                            mRData[outputPixelIdx] |= FLAG_OK_CONV;
                        }
                    }

                }
            }
        }
    }
#endif
    if (hasVarianceImage) {
        free(*variance);
        *variance = vData;
    }
    return;
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
double make_kernel(int xi, int yi, double *kernelSol) {

    int    gaussianCompIdx,solutionIdx,polyDegX,polyDegY,kernelPixelIdx;
    int    pixelCompIdx;
    double polyBasisX,polyBasisY,kernelSum;
    double normalizedX, normalizedY;

    solutionIdx  = 2;
    /* RANGE FROM -1 to 1 */
    normalizedX = (xi - 0.5 * rPixX) / (0.5 * rPixX);
    normalizedY = (yi - 0.5 * rPixY) / (0.5 * rPixY);

    for (gaussianCompIdx = 1; gaussianCompIdx < nCompKer; gaussianCompIdx++) {
        kernel_coeffs[gaussianCompIdx] = 0.0;
        polyBasisX = 1.0;
        for (polyDegX = 0; polyDegX <= kerOrder; polyDegX++) {
            polyBasisY = 1.0;
            for (polyDegY = 0; polyDegY <= kerOrder - polyDegX; polyDegY++) {
                kernel_coeffs[gaussianCompIdx] += kernelSol[solutionIdx++] * polyBasisX * polyBasisY;
                polyBasisY *= normalizedY;
            }
            polyBasisX *= normalizedX;
        }
    }
    kernel_coeffs[0] = kernelSol[1];

    for (kernelPixelIdx = 0; kernelPixelIdx < fwKernel * fwKernel; kernelPixelIdx++)
        kernel[kernelPixelIdx] = 0.0;

    kernelSum = 0.0;
    for (kernelPixelIdx = 0; kernelPixelIdx < fwKernel * fwKernel; kernelPixelIdx++) {
        for (pixelCompIdx = 0; pixelCompIdx < nCompKer; pixelCompIdx++) {
            kernel[kernelPixelIdx] += kernel_coeffs[pixelCompIdx] * kernel_vec[pixelCompIdx][kernelPixelIdx];
        }
        kernelSum += kernel[kernelPixelIdx];
    }
    return kernelSum;
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
double get_background(int xi, int yi, double *kernelSol) {

    double  background,polyBasisX,polyBasisY,normalizedX,normalizedY;
    int     polyDegX,polyDegY,solutionIdx;
    int     backgroundComponentOffset;

    backgroundComponentOffset = (nCompKer - 1) * ( ((kerOrder + 1) * (kerOrder + 2)) / 2 ) + 1;

    background = 0.0;
    solutionIdx     = 1;
    /* RANGE FROM -1 to 1 */
    normalizedX = (xi - 0.5 * rPixX) / (0.5 * rPixX);
    normalizedY = (yi - 0.5 * rPixY) / (0.5 * rPixY);

    polyBasisX=1.0;
    for (polyDegX = 0; polyDegX <= bgOrder; polyDegX++) {
        polyBasisY = 1.0;
        for (polyDegY = 0; polyDegY <= bgOrder - polyDegX; polyDegY++) {
            background += kernelSol[backgroundComponentOffset+solutionIdx++] * polyBasisX * polyBasisY;
            /* fprintf(stderr, "bg: %d %d %d %d %f %f %f\n", xi, yi, polyDegX, polyDegY, polyBasisX, polyBasisY, kernelSol[backgroundComponentOffset+solutionIdx-1]); */
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
void make_model(stamp_struct *stamp, double *kernelSol, float *csModel) {

    int       gaussianCompIdx,solutionIdx,polyDegX,polyDegY,stampPixelIdx,stampCenterX,stampCenterY;
    double    polyBasisX,polyBasisY,polynomialCoeff;
    double    *vector;
    float     halfPixX, halfPixY;
    double    normalizedX, normalizedY;

    halfPixX   = 0.5 * rPixX;
    halfPixY   = 0.5 * rPixY;

    stampCenterX = stamp->xss[stamp->sscnt];
    stampCenterY = stamp->yss[stamp->sscnt];

    /* RANGE FROM -1 to 1 */
    normalizedX = (stampCenterX - 0.5 * rPixX) / (0.5 * rPixX);
    normalizedY = (stampCenterY - 0.5 * rPixY) / (0.5 * rPixY);

    for (stampPixelIdx = 0; stampPixelIdx < fwKSStamp * fwKSStamp; stampPixelIdx++) csModel[stampPixelIdx] = 0.0;

    vector = stamp->vectors[0];
    polynomialCoeff  = kernelSol[1];
    for (stampPixelIdx = 0; stampPixelIdx < fwKSStamp * fwKSStamp; stampPixelIdx++) csModel[stampPixelIdx] += polynomialCoeff * vector[stampPixelIdx];

    solutionIdx=2;
    for (gaussianCompIdx = 1; gaussianCompIdx < nCompKer; gaussianCompIdx++) {
        vector = stamp->vectors[gaussianCompIdx];
        polynomialCoeff  = 0.0;
        polyBasisX     = 1.0;
        for (polyDegX = 0; polyDegX <= kerOrder; polyDegX++) {
            polyBasisY = 1.0;
            for (polyDegY = 0; polyDegY <= kerOrder - polyDegX; polyDegY++) {
                polynomialCoeff += kernelSol[solutionIdx++] * polyBasisX * polyBasisY;
                polyBasisY *= normalizedY;
            }
            polyBasisX *= normalizedX;
        }

        for (stampPixelIdx = 0; stampPixelIdx < fwKSStamp * fwKSStamp; stampPixelIdx++) {
            csModel[stampPixelIdx] += polynomialCoeff*vector[stampPixelIdx];
        }
    }
    return;
}

