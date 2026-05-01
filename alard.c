#include<stdio.h>
#include<string.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<fitsio.h>
#include<lapacke.h>

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
    int ig, idegx, idegy, nvec;
    int ren;
    
    nvec = 0;
    for (ig = 0; ig < ngauss; ig++) {
        for (idegx = 0; idegx <= deg_fixe[ig]; idegx++) {
            for (idegy = 0; idegy <= deg_fixe[ig]-idegx; idegy++) {
                /* stores kernel weight mask for each order */
                kernel_vec[nvec] = kernel_vector(nvec, idegx, idegy, ig, &ren);
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
    
    int       ren = 0;
    int       i,j,xi,yi,dx,dy,idegx,idegy,di,dj,nv,ig,nvec;
    double    ax,ay,xf,yf;
    double *im;
    float     rPixX2, rPixY2;
    
    rPixX2   = 0.5 * rPixX;
    rPixY2   = 0.5 * rPixY;
    
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
    
    nvec = 0;
    for (ig = 0; ig < ngauss; ig++) {
        for (idegx = 0; idegx <= deg_fixe[ig]; idegx++) {
            for (idegy = 0; idegy <= deg_fixe[ig]-idegx; idegy++) {
                
                ren = 0;
                dx = (idegx / 2) * 2 - idegx;
                dy = (idegy / 2) * 2 - idegy;
                if (dx == 0 && dy == 0 && nvec > 0)
                    ren = 1;
                
                /* fill stamp->vectors[nvec] with convolved image */
                /* image is convolved with functional form of kernel, fit later for amplitude */
                xy_conv_stamp(stamp, imConv, nvec, ren);
                ++nvec;
            }
        }
    }
    
    /* get the krefArea data */
    if (cutSStamp(stamp, imRef))
        return 1;
    
    /* fill stamp->vectors[nvec+++] with x^(bg) * y^(bg) for background fit */
    xi = stamp->xss[stamp->sscnt];
    yi = stamp->yss[stamp->sscnt];
    di = xi - hwKSStamp;
    dj = yi - hwKSStamp;
    for (i = xi - hwKSStamp; i <= xi + hwKSStamp; i++) {
        xf = (i - rPixX2) / rPixX2;
        
        for (j = yi - hwKSStamp; j <= yi + hwKSStamp; j++) {
            /* fprintf(stderr, "%d %d %d %d %d %d\n", k, xi, yi,i, j, fwKSStamp); */
            yf = (j - rPixY2) / rPixY2; 
            
            ax = 1.0;
            nv = nvec;
            for (idegx = 0; idegx <= bgOrder; idegx++) {
                ay = 1.0; 
                for (idegy = 0; idegy <= bgOrder - idegx; idegy++) {
                    im = stamp->vectors[nv];
                    im[i-di+fwKSStamp*(j-dj)] = ax * ay;
                    ay *= yf;
                    ++nv;
                }
                ax *= xf;
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
    int       i,j,k,dx,dy,ix;
    double    sum_x,sum_y,x,qe;
    
    if (usePCA) {
        return kernel_vector_PCA(n, deg_x, deg_y, ig, ren);
    }
    
    vector = (double *)malloc(fwKernel*fwKernel*sizeof(double));
    dx = (deg_x / 2) * 2 - deg_x;
    dy = (deg_y / 2) * 2 - deg_y;
    sum_x = sum_y = 0.0;
    *ren = 0;
    
    for (ix = 0; ix < fwKernel; ix++) {
        x            = (double)(ix - hwKernel);
        k            = ix+n*fwKernel;
        qe           = exp(-x * x * sigma_gauss[ig]);
        filter_x[k]  = qe * pow(x, deg_x);
        filter_y[k]  = qe * pow(x, deg_y);
        sum_x       += filter_x[k];
        sum_y       += filter_y[k];
    }
    
    if (n > 0)
        kernel0 = kernel_vec[0];
    
    sum_x = 1. / sum_x;
    sum_y = 1. / sum_y;
    
    if (dx == 0 && dy == 0) {
        for (ix = 0; ix < fwKernel; ix++) {
            filter_x[ix+n*fwKernel] *= sum_x;
            filter_y[ix+n*fwKernel] *= sum_y;
        }
        
        for (i = 0; i < fwKernel; i++) {
            for (j = 0; j < fwKernel; j++) {
                vector[i+fwKernel*j] = filter_x[i+n*fwKernel] * filter_y[j+n*fwKernel];
            }
        }
        
        if (n > 0) {
            for (i = 0; i < fwKernel * fwKernel; i++) {
                vector[i] -= kernel0[i];
            }
            *ren = 1;
        }
    } else {
        for (i = 0; i < fwKernel; i++) {
            for (j = 0; j < fwKernel; j++) {
                vector[i+fwKernel*j] = filter_x[i+n*fwKernel] * filter_y[j+n*fwKernel];
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
    int i,j;
    
    vector = (double *)malloc(fwKernel*fwKernel*sizeof(double));
    
    for (i = 0; i < fwKernel; i++) {
        for (j = 0; j < fwKernel; j++) {
            vector[i+fwKernel*j] = PCA[n][i+fwKernel*j];
        }
    }
    
    if (n > 0)
        kernel0 = kernel_vec[0];
    
    if (n > 0) {
        for (i = 0; i < fwKernel * fwKernel; i++) {
            vector[i] -= kernel0[i];
        }
        *ren = 1;
    }
    
    return vector;
}

/**
 * @brief Convolve the substamp region of an image with the n-th separable
 *        Gaussian-polynomial kernel basis function.
 *
 * @details Uses the separability of the Gaussian-polynomial filter to perform
 * the convolution in two 1-D passes (first along y, then along x), storing the
 * fwKSStamp×fwKSStamp result in stamp->vectors[n].  This is the performance-
 * critical inner loop (~60 % of total runtime according to the PROF file).
 *
 * If the renormalise flag @p ren is set (i.e. kernel_vector() returned ren=1
 * for this basis), the zeroth-basis convolved image stamp->vectors[0] is
 * subtracted pixel-by-pixel from the result, matching the differential
 * construction of the higher-order basis images.
 *
 * In PCA mode the call is forwarded to xy_conv_stamp_PCA().
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
    int       i,j,xc,yc,xij,sub_width,xi,yi;
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
    for(i = xi - hwKSStamp - hwKernel; i <= xi + hwKSStamp + hwKernel; i++) {
        for(j = yi - hwKSStamp; j <= yi + hwKSStamp; j++) {
            xij = i - xi + sub_width / 2 + sub_width * (j - yi + hwKSStamp);
            temp[xij] = 0.0;
            for(yc = -hwKernel; yc <= hwKernel; yc++) {
                temp[xij] += image[i+rPixX*(j+yc)] * filter_y[hwKernel-yc+n*fwKernel];
            }
        }
    }
    
    /* convolve with x filter */
    for(j = -hwKSStamp; j <= hwKSStamp; j++) {
        for(i = -hwKSStamp; i <= hwKSStamp;i++) {  
            xij = i + hwKSStamp + fwKSStamp * (j + hwKSStamp);
            imc[xij] = 0.0;
            for(xc = -hwKernel; xc <= hwKernel; xc++) {
                imc[xij] += temp[i+xc+sub_width/2+sub_width*(j+hwKSStamp)] * filter_x[hwKernel-xc+n*fwKernel];
            }
        }
    }
    
    if (ren) {
        v0 = stamp->vectors[0];
        for(i = 0; i < fwKSStamp * fwKSStamp; i++) imc[i] -= v0[i];
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
    
    int       i,j,xc,yc,xij,xi,yi;
    double    *v0,*imc;
    
    xi  = stamp->xss[stamp->sscnt];
    yi  = stamp->yss[stamp->sscnt];
    imc = stamp->vectors[n];
    
    /* pull area to convolve out of full reference image region */
    for(j = yi - hwKSStamp; j <= yi + hwKSStamp; j++) {
        for(i = xi - hwKSStamp; i <= xi + hwKSStamp; i++) {
            xij      = i - (xi - hwKSStamp) + fwKSStamp * (j - (yi - hwKSStamp));
            imc[xij] = 0.;
            
            for(yc = -hwKernel; yc <= hwKernel; yc++) {
                for(xc = -hwKernel; xc <= hwKernel; xc++) {
                    imc[xij] += image[(i+xc)+rPixX*(j+yc)] * PCA[n][(xc+hwKernel) + fwKernel*(yc+hwKernel)];
                }
            }
        }
    }
    
    if (ren) {
        v0 = stamp->vectors[0];
        for(i = 0; i < fwKSStamp * fwKSStamp; i++) imc[i] -= v0[i];
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
    char check;
    int i,mat_size;
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
    for (i = 0; i <= mat_size; i++) 
        matrix[i] = (double *)malloc((mat_size + 1)*sizeof(double));
    
    /* allocate weight matrix */
    wxy = (double **)malloc(nS*sizeof(double *));
    for (i = 0; i < nS; i++)
        wxy[i] = (double *)malloc(ncomp2*sizeof(double));
    
    
    
    if (verbose>=2) fprintf(stderr, " Expanding Matrix For Full Fit\n");
    build_matrix(stamps, nS, matrix);
    build_scprod(stamps, nS, imRef, kernelSol);
    
    solve_spd(matrix, mat_size, kernelSol);

    if (verbose>=2) fprintf(stderr, " Checking again\n");
    check = check_again(stamps, kernelSol, imConv, imRef, imNoise, meansigSubstamps, scatterSubstamps, NskippedSubstamps);

    while(check) {

        fprintf(stderr, "\n Re-Expanding Matrix\n");
        build_matrix(stamps, nS, matrix);
        build_scprod(stamps, nS, imRef, kernelSol);

        solve_spd(matrix, mat_size, kernelSol);
        
        fprintf(stderr, " Checking again\n");          
        check = check_again(stamps, kernelSol, imConv, imRef, imNoise, meansigSubstamps, scatterSubstamps, NskippedSubstamps); 
    }
    fprintf(stderr, " Sigma clipping of bad stamps converged, kernel determined\n");
    
    for (i = 0; i <= mat_size; i++)
        free(matrix[i]);
    for (i = 0; i < nS; i++)
        free(wxy[i]);
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
    
    int       i,j,pixStamp,k,i1,ivecbg=0;
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
    for (i = 0; i < ncomp1; i++) {
        for (j = 0; j <= i; j++) {
            q = 0.0;
            /* integrate W_m1 and W_m2 (sum over all pixels) */
            for (k = 0; k < pixStamp; k++) 
                q += vec[i][k] * vec[j][k];
            
            /* Q from Eqn 3. in Alard */
            stamp->mat[i+1][j+1] = q;
        }
    }
    
    for (i1 = 0; i1 < ncomp1; i1++) { 
        ivecbg = ncomp1;
        
        p0 = 0.0;
        /* integrate convolved images and first order background (equals 1 everywhere!)*/
        for (k = 0; k < pixStamp; k++)
            p0 += vec[i1][k] * vec[ivecbg][k];
        stamp->mat[ncomp1+1][i1+1] = p0;
    }  
    
    /* integrate first order background with itself */
    /* NOTE : DON'T MASK K HERE - BACKGROUND! */
    for (k = 0, q = 0.0; k < pixStamp; k++)
        q += vec[ivecbg][k] * vec[ncomp1][k];
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
    
    int       xc,yc,xi,yi,i1,k;
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
    for (i1 = 0; i1 < ncomp1; i1++) {
        p0 = 0.0;
        for (xc = -hwKSStamp; xc <= hwKSStamp; xc++) {
            for (yc = -hwKSStamp; yc <= hwKSStamp; yc++) {
                k   = xc + hwKSStamp + fwKSStamp * (yc + hwKSStamp);
                p0 += vec[i1][k] * image[xc+xi+rPixX*(yc+yi)];
            }
        }
        stamp->scprod[i1+1] = p0;
    }
    
    /* Multiply first order background model with reference image */
    q = 0.0;
    for (xc = -hwKSStamp; xc <= hwKSStamp; xc++) {
        for (yc = -hwKSStamp; yc <= hwKSStamp; yc++) {
            k  = xc + hwKSStamp + fwKSStamp * (yc + hwKSStamp);
            q += vec[ncomp1][k] * image[xc+xi+rPixX*(yc+yi)];	 
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
    
    int    nComps,i,im,jm,mcnt1,mcnt2,mcnt3;
    double sum=0,kmean,kstdev;
    double merit1,merit2,merit3,sig1,sig2,sig3;
    float *m1,*m2,*m3,*ks;
    int    xc, yc, nks;
    
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
    
    for (i = 0; i < nS; i++) {
        
        xc    = stamps[i].xss[stamps[i].sscnt];
        yc    = stamps[i].yss[stamps[i].sscnt];
        
        /* extract check_mat to solve one particular stamp */
        for (im = 1; im <= nComps; im++) {
            check_vec[im] = stamps[i].scprod[im];
            
            for (jm = 1; jm <= im; jm++) {
                check_mat[im][jm] = stamps[i].mat[im][jm];
                check_mat[jm][im] = check_mat[im][jm];
            }
        }
        
        /* fit stamp, the constant kernel coefficients end up in check_vec */
        solve_spd(check_mat, nComps, check_vec);
        
        /* find kernel sum */
        sum = check_vec[1];
        check_stack[i] = sum;
        stamps[i].norm = sum;
        ks[nks++]      = sum;
        
        if (verbose >= 2) fprintf(stderr, "    # %d    xss: %4i yss: %4i  ksum: %f\n", i,
                                  stamps[i].xss[stamps[i].sscnt],
                                  stamps[i].yss[stamps[i].sscnt], sum);
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
    for (i = 0; i < nS; i++) {
        stamps[i].diff = fabs((stamps[i].norm - kmean) / kstdev);
    }
    
    /*****************************************************
     * Global fit for kernel solution
     *****************************************************/
    
    /* do only if necessary */
    if ((strncmp(forceConvolve, "b", 1)==0)) {
        
        /* allocate fitting matrix */
        matrix = (double **)calloc((mat_size + 1), sizeof(double *));
        for (i = 0; i <= mat_size; i++) 
            matrix[i] = (double *)calloc((mat_size + 1), sizeof(double));
        
        /* allocate weight matrix */
        wxy = (double **)calloc(nS, sizeof(double *));
        for (i = 0; i < nS; i++)
            wxy[i] = (double *)calloc(ncomp2, sizeof(double));
        
        /* first find out how many good stamps to allocate */
        ntestStamps = 0;
        for (i = 0; i < nS; i++)
            if (stamps[i].diff < kerSigReject) {
                ntestStamps++;
            }
            else {
                if (verbose >= 2) fprintf(stderr, "    # %d    skipping xss: %4i yss: %4i ksum: %f sigma: %f\n", i,
                                          stamps[i].xss[stamps[i].sscnt],
                                          stamps[i].yss[stamps[i].sscnt],
                                          stamps[i].norm, stamps[i].diff);
            }
        
        /* then allocate test stamp structure */
        if(!(testStamps = (stamp_struct *)calloc(ntestStamps, sizeof(stamp_struct)))) {
            printf("Cannot Allocate Test Stamp List\n"); 
            exit (1);
        }
        testKerSol = (double *)calloc((nCompTotal+1), sizeof(double));
        
        /* and point test stamp structure to good stamps */
        ntestStamps = 0;
        for (i = 0; i < nS; i++)
            if (stamps[i].diff < kerSigReject)
                testStamps[ntestStamps++] = stamps[i];
        
        /* finally do fit */
        if (verbose >= 2) fprintf(stderr, " Expanding Test Matrix For Fit\n");
        build_matrix(testStamps, ntestStamps, matrix);
        build_scprod(testStamps, ntestStamps, imRef, testKerSol);
        solve_spd(matrix, mat_size, testKerSol);
        
        /* get the kernel sum to normalize figures of merit! */
        kmean = make_kernel(0, 0, testKerSol);
        
        /* determine figure of merit from good stamps */
        
        /* average of sum (diff**2 / value), ~variance */
        m1 = (float *)calloc(ntestStamps, sizeof(float));
        
        /* standard deviation of pixel distribution */
        m2 = (float *)calloc(ntestStamps, sizeof(float));
        
        /* noise sd based on histogram distribution width */
        m3 = (float *)calloc(ntestStamps, sizeof(float));
        
        mcnt1 = 0;
        mcnt2 = 0;
        mcnt3 = 0;
        for (i = 0; i < ntestStamps; i++) {
            
            getStampSig(&testStamps[i], testKerSol, imNoise, &sig1, &sig2, &sig3);
            
            if ((sig1 != -1) && (sig1 <= MAXVAL)) {
                m1[mcnt1++] = sig1;
            }
            if ((sig2 != -1) && (sig2 <= MAXVAL)) {
                m2[mcnt2++] = sig2;
            }
            if ((sig3 != -1) && (sig3 <= MAXVAL)) {
                m3[mcnt3++] = sig3;
            }
        }
        sigma_clip(m1, mcnt1, &merit1, &sig1, 10);
        sigma_clip(m2, mcnt2, &merit2, &sig2, 10);
        sigma_clip(m3, mcnt3, &merit3, &sig3, 10);
        
        /* normalize by kernel sum */
        merit1 /= kmean;
        merit2 /= kmean;
        merit3 /= kmean;
        
        /* clean up this mess */
        if (testKerSol)	 free(testKerSol);
        if (testStamps)	 free(testStamps);
        for (i = 0; i <= mat_size; i++)
            free(matrix[i]);
        for (i = 0; i < nS; i++)
            free(wxy[i]);
        free(matrix);
        free(wxy);
        
        free(m1);
        free(m2);
        free(m3);
        free(ks);
        
        /* average value of figures of merit across stamps */
        fprintf(stderr, "    <var_merit> = %.3f, <sd_merit> = %.3f, <hist_merit> = %.3f\n", merit1, merit2, merit3);
        
        /* return what is asked for if possible, if not use backup */
        if (strncmp(figMerit, "v", 1)==0) {
            if (mcnt1 > 0) {
                return merit1;
            }
            else if (mcnt2 > 0) {
                return merit2;
            }
            else if (mcnt3 > 0) {
                return merit3;
            }
            else {
                return 666;
            }
        }
        else if (strncmp(figMerit, "s", 1)==0) {
            if (mcnt2 > 0) {
                return merit2;
            }
            else if (mcnt1 > 0) {
                return merit1;
            }
            else if (mcnt3 > 0) {
                return merit3;
            }
            else {
                return 666;
            }
        }
        else if (strncmp(figMerit, "h", 1)==0) {      
            if (mcnt3 > 0) {
                return merit3;
            }
            else if (mcnt1 > 0) {
                return merit1;
            }
            else if (mcnt2 > 0) {
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
            i1 = i / ncomp2;
            i2 = i - i1 * ncomp2;
            
            for (j = 0; j <= i; j++) {
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
    
    int       mat_size,i,j,pixStamp,istamp,k,i1,i2,j1,j2,ibg,jbg,ivecbg,jj;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    **matrix0,p0,q;
    double    **vec;
    float     rPixX2, rPixY2;
    
    int ideg1, ideg2, xstamp, ystamp;
    double a1, a2, fx, fy;
    
    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    pixStamp = fwKSStamp * fwKSStamp;
    rPixX2   = 0.5 * rPixX;
    rPixY2   = 0.5 * rPixY;
    
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
        /* RANGE FROM -1 to 1 */
        fx     = (xstamp - rPixX2) / rPixX2; 
        fy     = (ystamp - rPixY2) / rPixY2; 
        
        /* build weight function *HERE* */
        k  = 0;
        a1 = 1.0;
        for (ideg1 = 0; ideg1 <= kerOrder; ideg1++) {
            a2 = 1.0;
            for (ideg2 = 0; ideg2 <= kerOrder - ideg1; ideg2++) {
                wxy[istamp][k++] = a1 * a2;
                a2 *= fy;
            }
            a1 *= fx;
        }
        
        matrix0 = stamps[istamp].mat;
        for (i = 0; i < ncomp; i++) {
            i1 = i / ncomp2;
            i2 = i - i1 * ncomp2;
            
            for (j = 0; j <= i; j++) {
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
    
    int       istamp,xc,yc,xi,yi,i1,i2,k,ibg,i,ii;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    p0,q;
    double    **vec;
    
    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    
    for (i = 0; i <= ncomp + nbg_vec + 1; i++)
        kernelSol[i]=0.0;
    
    for (istamp = 0; istamp < nS; istamp++) {
        /* skip over any bad stamps along the way */
        while(stamps[istamp].sscnt >= stamps[istamp].nss) {
            ++istamp; 
            if (istamp >= nS) break;
        }
        if (istamp >= nS) break;
        
        vec= stamps[istamp].vectors;
        xi = stamps[istamp].xss[stamps[istamp].sscnt];
        yi = stamps[istamp].yss[stamps[istamp].sscnt];
        
        p0 = stamps[istamp].scprod[1];
        kernelSol[1] += p0;
        
        /* spatially weighted convolved image * ref image */
        for (i1 = 1; i1 < ncomp1 + 1; i1++) {
            p0 = stamps[istamp].scprod[i1+1];
            for (i2 = 0; i2 < ncomp2; i2++) {
                ii = (i1-1) * ncomp2 + i2 + 1;
                /* no need for weighting here */
                kernelSol[ii+1] += p0 * wxy[istamp][i2];
            }
        }
        
        /* spatially weighted bg model convolved with ref image */
        for (ibg = 0; ibg < nbg_vec; ibg++) {
            q = 0.0;
            for (xc = -hwKSStamp; xc <= hwKSStamp; xc++) {
                for (yc = -hwKSStamp; yc <= hwKSStamp;yc++) {
                    k  = xc + hwKSStamp + fwKSStamp * (yc + hwKSStamp);
                    q += vec[ncomp1+ibg+1][k] * image[xc+xi+rPixX*(yc+yi)];
                    
                }
            }
            kernelSol[ncomp+ibg+2] += q;
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
    int i, j, idx, nsig=0;
    int xRegion2, xRegion = stamp->xss[stamp->sscnt];
    int yRegion2, yRegion = stamp->yss[stamp->sscnt];
    float idat, indat;
    
    *sig = 0;
    
    for (j = 0; j < fwKSStamp; j++) {
        yRegion2 = yRegion - hwKSStamp + j;
        
        for (i = 0; i < fwKSStamp; i++) {
            xRegion2 = xRegion - hwKSStamp + i;
            
            idx   = xRegion2+rPixX*yRegion2;
            idat  = imDiff[idx];
            indat = 1. / imNoise[idx];
            
            /* this shouldn't be the case, but just in case... */
            if (mRData[idx] & FLAG_INPUT_ISBAD)
                continue;
            
            nsig++;
            *sig += idat * idat * indat * indat;
            
        }
    }
    if (nsig > 0) 
        *sig /= nsig;
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
    int i, j, idx, sscnt, nsig, xRegion, yRegion, xRegion2, yRegion2;
    double cSum, cMean, cMedian, cMode, cLfwhm;
    double *im, tdat, idat, ndat, diff, bg;
    
    /* info */
    sscnt   = stamp->sscnt;
    xRegion = stamp->xss[sscnt];
    yRegion = stamp->yss[sscnt];
    
    /* the comparison image */
    im = stamp->krefArea;
    /* background from fit */
    bg = get_background(xRegion, yRegion, kernelSol);
    /* temp contains the convolved image from fit, fwKSStamp x fwKSStamp */
    make_model(stamp, kernelSol, temp); 
    
    /* get sigma of stamp diff */
    nsig = 0;
    *sig1 = 0;
    *sig2 = 0;
    *sig3 = 0;
    
    for (j = 0; j < fwKSStamp; j++) {
        yRegion2 = yRegion - hwKSStamp + j;
        
        for (i = 0; i < fwKSStamp; i++) {
            xRegion2 = xRegion - hwKSStamp + i;
            
            idx  = i+j*fwKSStamp;
            
            tdat = temp[idx];
            idat = im[idx];
            ndat = imNoise[xRegion2+rPixX*yRegion2];
            
            diff = tdat - idat + bg;
            
            if ((mRData[xRegion2+rPixX*yRegion2] & FLAG_INPUT_ISBAD) || (fabs(idat) <= ZEROVAL)) {
                continue;
            }
            else {
                temp[idx] = diff;
            }
            
            /* check for NaN */
            if ((tdat*0.0 != 0.0) || (idat*0.0 != 0.0)) {
                mRData[xRegion2+rPixX*yRegion2] |= (FLAG_INPUT_ISBAD | FLAG_ISNAN);
                continue;
            }
            
            nsig++;
            *sig1 += diff * diff / ndat;
            /*fprintf(stderr, "OK %d %d : %f %f %f\n", xRegion2, yRegion2, tdat, idat, ndat);*/
        }
    }
    if (nsig > 0) {
        *sig1 /= nsig;
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
void spatial_convolve(float *image, float **variance, int xSize, int ySize, double *kernelSol, float *cRdata, int *cMask) {
    
    int       i1,j1,i2,j2,nsteps_x,nsteps_y,i,j,i0,j0,ic,jc,ik,jk,nc,ni,mbit,dovar;
    double    q, qv, kk, aks, uks;
    float     *vData=NULL;
    
    if ((*variance) == NULL)
        dovar = 0;
    else
        dovar = 1;
    
    if (dovar) {
        if ( !(vData = (float *)calloc(xSize*ySize, sizeof(float)))) {
            return;
        }
    }
    
    nsteps_x = ceil((double)(xSize)/(double)kcStep);
    nsteps_y = ceil((double)(ySize)/(double)kcStep);
    
    for (j1 = 0; j1 < nsteps_y; j1++) {
        j0 = j1 * kcStep + hwKernel;
        
        for(i1 = 0; i1 < nsteps_x; i1++) {
            i0 = i1 * kcStep + hwKernel;
            
            make_kernel(i0 + hwKernel, j0 + hwKernel, kernelSol);
            
            for (j2 = 0; j2 < kcStep; j2++) {
                j = j0 + j2;
                if ( j >= ySize - hwKernel) break;
                
                for (i2 = 0; i2 < kcStep; i2++) {
                    i = i0 + i2;
                    if (i >= xSize - hwKernel) break;
                    
                    ni = i+xSize*j;
                    qv = q = aks = uks = 0.0;
                    mbit = 0x0;
                    for (jc = j - hwKernel; jc <= j + hwKernel; jc++) {
                        jk = j - jc + hwKernel;
                        
                        for (ic = i - hwKernel; ic <= i + hwKernel; ic++) {
                            ik = i - ic + hwKernel;
                            
                            nc = ic+xSize*jc;
                            kk = kernel[ik+jk*fwKernel];
                            
                            q     += image[nc] * kk;
                            if (dovar) {
                                if (convolveVariance)
                                    qv += (*variance)[nc] * kk;
                                else
                                    qv += (*variance)[nc] * kk * kk;			   
                            }
                            
                            mbit  |= cMask[nc];
                            aks   += fabs(kk);
                            if (!(cMask[nc] & FLAG_INPUT_ISBAD)) {
                                uks += fabs(kk);
                            }
                        }
                    }
                    
                    cRdata[ni]   = q;
                    if (dovar)
                        vData[ni] = qv;
                    
                    /* mask propagation changed in 5.1.9 */
                    /* mRData[ni]  |= mbit; */ 
                    /* mRData[ni]  |= FLAG_OK_CONV      * (mbit > 0);*/
                    mRData[ni]  |= cMask[ni];
                    mRData[ni]  |= FLAG_OUTPUT_ISBAD * ((cMask[ni] & FLAG_INPUT_ISBAD) > 0);
                    
                    if (mbit) {
                        if ((uks / aks) < kerFracMask) {
                            mRData[ni] |= (FLAG_OUTPUT_ISBAD | FLAG_BAD_CONV);
                        }
                        else {
                            mRData[ni] |= FLAG_OK_CONV;
                        }
                    }
                    
                }       
            }
        }
    }
    if (dovar) {
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
    
    int    i1,k,ix,iy,i;
    double ax,ay,sum_kernel;
    double xf, yf;
    
    k  = 2;
    /* RANGE FROM -1 to 1 */
    xf = (xi - 0.5 * rPixX) / (0.5 * rPixX);
    yf = (yi - 0.5 * rPixY) / (0.5 * rPixY);
    
    for (i1 = 1; i1 < nCompKer; i1++) {
        kernel_coeffs[i1] = 0.0;
        ax = 1.0;
        for (ix = 0; ix <= kerOrder; ix++) {
            ay = 1.0;
            for (iy = 0; iy <= kerOrder - ix; iy++) {
                kernel_coeffs[i1] += kernelSol[k++] * ax * ay;
                ay *= yf;
            }
            ax *= xf;
        }
    }
    kernel_coeffs[0] = kernelSol[1]; 
    
    for (i = 0; i < fwKernel * fwKernel; i++)
        kernel[i] = 0.0;
    
    sum_kernel = 0.0;
    for (i = 0; i < fwKernel * fwKernel; i++) {
        for (i1 = 0; i1 < nCompKer; i1++) {
            kernel[i] += kernel_coeffs[i1] * kernel_vec[i1][i];
        }
        sum_kernel += kernel[i];    
    }
    return sum_kernel;
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
    
    double  background,ax,ay,xf,yf;
    int     i,j,k;
    int     ncompBG;
    
    ncompBG = (nCompKer - 1) * ( ((kerOrder + 1) * (kerOrder + 2)) / 2 ) + 1;
    
    background = 0.0;
    k          = 1;
    /* RANGE FROM -1 to 1 */
    xf = (xi - 0.5 * rPixX) / (0.5 * rPixX);
    yf = (yi - 0.5 * rPixY) / (0.5 * rPixY);
    
    ax=1.0;
    for (i = 0; i <= bgOrder; i++) {
        ay = 1.0; 
        for (j = 0; j <= bgOrder - i; j++) {
            background += kernelSol[ncompBG+k++] * ax * ay;
            /* fprintf(stderr, "bg: %d %d %d %d %f %f %f\n", xi, yi, i, j, ax, ay, kernelSol[ncompBG+k-1]); */
            ay *= yf;
        }
        ax *= xf;
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
    
    int       i1,k,ix,iy,i,xi,yi;
    double    ax,ay,coeff;
    double    *vector;
    float     rPixX2, rPixY2;
    double    xf, yf;
    
    rPixX2   = 0.5 * rPixX;
    rPixY2   = 0.5 * rPixY;
    
    xi = stamp->xss[stamp->sscnt];
    yi = stamp->yss[stamp->sscnt];
    
    /* RANGE FROM -1 to 1 */
    xf = (xi - 0.5 * rPixX) / (0.5 * rPixX);
    yf = (yi - 0.5 * rPixY) / (0.5 * rPixY);
    
    for (i = 0; i < fwKSStamp * fwKSStamp; i++) csModel[i] = 0.0;
    
    vector = stamp->vectors[0];
    coeff  = kernelSol[1];
    for (i = 0; i < fwKSStamp * fwKSStamp; i++) csModel[i] += coeff * vector[i];
    
    k=2;
    for (i1 = 1; i1 < nCompKer; i1++) {
        vector = stamp->vectors[i1];
        coeff  = 0.0; 
        ax     = 1.0;
        for (ix = 0; ix <= kerOrder; ix++) {
            ay = 1.0;
            for (iy = 0; iy <= kerOrder - ix; iy++) {
                coeff += kernelSol[k++] * ax * ay;
                ay *= yf;
            }
            ax *= xf;
        }
        
        for (i = 0; i < fwKSStamp * fwKSStamp; i++) {
            csModel[i] += coeff*vector[i];
        }
    }
    return;
}

