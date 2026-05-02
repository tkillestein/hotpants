#define FLAG_BAD_PIXVAL         0x01   /* 1 */
#define FLAG_SAT_PIXEL          0x02   /* 2 */
#define FLAG_LOW_PIXEL          0x04   /* 4 */
#define FLAG_ISNAN              0x08   /* 8 */

#define FLAG_BAD_CONV           0x10   /* 16  */
#define FLAG_INPUT_MASK         0x20   /* 32  */
#define FLAG_OK_CONV            0x40   /* 64  */
#define FLAG_INPUT_ISBAD        0x80   /* 128 */

#define FLAG_T_BAD              0x100  /* 256  */
#define FLAG_T_SKIP             0x200  /* 512  */
#define FLAG_I_BAD              0x400  /* 1024 */
#define FLAG_I_SKIP             0x800  /* 2048 */

#define FLAG_OUTPUT_ISBAD       0x8000 /* 32768 */


#define MAXDIM        4
#define SCRLEN        256
#define MAXVAL        1e10
#define ZEROVAL       1e-10

#define D_NGAUSS      3       
#define D_DEG_GAUSS1  6       
#define D_DEG_GAUSS2  4       
#define D_DEG_GAUSS3  2       
#define D_SIG_GAUSS1  0.7     
#define D_SIG_GAUSS2  1.5     
#define D_SIG_GAUSS3  3.0

#define D_UTHRESH       25000.  /* assumed upper valid data */
#define D_LTHRESH       0.      /* assumed lower valid data */
#define D_GAIN          1.      /* assumed gain, e-/ADU */
#define D_RDNOISE       0.      /* assumed readnoise, e- */
#define D_PEDESTAL      0.      /* assumed pedestal in ADU */
#define D_HWKERNEL      10      /* pixel half width of kernel */

#define D_KFITTHRESH    20.     /* # sigma above noise to use in kernel fit */
#define D_NFITTHRESH    0.1     /* if this fraction of requested stamps are not obtained... */
#define D_SFITTHRESH    0.5     /* scale D_KFITTHRESH down by this factor */

#define D_INMASKFSPREAD 1.      /* fraction of kernel width to spread input mask */
#define D_OUMASKFSPREAD 1.      /* fraction of kernel width to spread output mask */
#define D_NREGIONS      1       /* number of regions per image */
#define D_NSTAMPS       10      /* number of stamps each region, along each axis */
#define D_USEFULLSS     0       /* maximally divide the image up into stamps, each the size of kernel */
#define D_FINDSSCENTERS 1       /* if file submitted containing substamps, fill rest of centers automatically? */
#define D_CMPFILE       0       /* if ssfile, is it a cmp file? */
#define D_NKSSTAMPS     3       /* number of kernel test substamps each stamp */
#define D_HWKSSTAMP     15      /* half size of kernel substamp in each stamp */
#define D_KORDER        2       /* order of spatially varying kernel */
#define D_BGORDER       1       /* order of spatially varying sky */
#define D_KSIGREJECT    2.      /* reject stamp at this sigma from kernel fit */
#define D_STATSIG       3.      /* threshold in sigma clipping algorithms */
#define D_KFRACMASK     0.99    /* fraction of abs(kernel) values below which mask bit 0x8000 is set */
#define D_RESCALEOK     0       /* rescale noise of 'OK' pixels */
#define D_NORMALIZE     "t"     /* output image on template phot system */
#define D_CONVOLVE      "b"     /* try both ways by default */
#define D_FIGMERIT      "v"     /* (v)ariance, (s)igma based or (h)istogram based method to direct convolution */
#define D_FILL          1.e-30  /* value put on masked pixels */
#define D_FILLNOISE     0       /* value put on masked pixels in noise image */
#define D_NOILAYER      0       /* add pure noise image as output image layer */
#define D_SIGLAYER      0       /* add noise-scaled diff image as output image layer */
#define D_CONVLAYER     0       /* add convolved image as output image layer */
#define D_OUTSUM        0       /* output image is *SUM* of input images */
#define D_NOCLOBBER     0       /* do not clobber output image - nice double negative */
#define D_KINFO         0       /* write out kernel info to output diff image - lots-o-keywords */
#define D_OUTSHORT      0       /* write out images as short ints (16 bitpix) instead of -32 bitpix */
#define D_OUTBZERO      0.      /* if OUTSHORT, use this bzero for output images */
#define D_OUTBSCALE     1.      /* if OUTSHORT, use this bscale for output images */
#define D_SAVEXY        0       /* armin stuff : save xy location of stamps */
#define D_VERBOSE       1       /* levels of verbosity, 0, 1, 2 */
#define D_CONVVAR       0       /* instead of convolving noise, convolve variance.  ker vs ker**2 */
#define D_NTHREAD       1       /* number of threads to use for calculating regions - will be limited to nrx*nry */

/* =====================================================================
   HISTOGRAM STATISTICS PARAMETERS (used in getStampStats3)
   =====================================================================
   getStampStats3 estimates image statistics (mean, median, mode, FWHM)
   by building a 256-bin histogram and performing sigma-clipping.
   See functions.c:getStampStats3() for algorithm details.
   ===================================================================== */

#define HISTOGRAM_SAMPLE_SIZE    100    /* # random pixels sampled to estimate bin width; must be << total pixels */
#define HISTOGRAM_NUM_BINS       256    /* # histogram bins for FWHM estimation (256 fits percentile calculations) */
#define HISTOGRAM_LOWER_FRAC     0.5    /* lower percentile (50%) for bin-width estimation */
#define HISTOGRAM_UPPER_FRAC     0.9    /* upper percentile (90%) for bin-width estimation */
#define HISTOGRAM_PEAK_WIDTH_PCT 0.1    /* width of histogram peak containing ~10% of points (used for mode finding) */
#define HISTOGRAM_NOISE_HALF_PCT 0.25   /* half-width percentile for FWHM: finds 25th-75th percentile range around mode */

/* =====================================================================
   RANDOM NUMBER GENERATION
   ===================================================================== */
#define RNG_SEED_MAGIC           (-666)  /* Numerical Recipes convention; Park & Miller generator re-seeds on negative value */

/* =====================================================================
   KERNEL-FITTING ITERATION LIMITS
   ===================================================================== */
#define MAX_SIGMA_CLIP_RETRIES   5       /* max iterations for histogram bin-width refinement in getStampStats3; typical convergence < 3 */

/* =====================================================================
   POLYNOMIAL BASIS FORMULA
   =====================================================================
   The number of polynomial basis terms (1-D and 2-D) up to order n is:
     N(n) = (n+1)*(n+2)/2

   This formula counts all monomials x^i*y^j where i,j >= 0 and i+j <= n.
   Used to compute:
   - nComp = number of spatial polynomial terms in kernel
   - ncompBG = number of background polynomial terms
   - nBGVectors = number of background basis vectors

   Reference: Alard & Lupton (1998), Eq. 3; Appendix A discusses the
   polynomial expansion of the spatially-varying kernel coefficient c_i(x,y).

   Example: order=2 gives (2+1)*(2+2)/2 = 6 terms: {1, x, y, x^2, xy, y^2}
   ===================================================================== */

/* Compile-time guards for type-size assumptions (C11/C17 _Static_assert) */
_Static_assert(sizeof(float)  == 4, "CFITSIO BITPIX=-32 I/O assumes 4-byte float");
_Static_assert(sizeof(double) == 8, "LAPACK/BLAS routines assume 8-byte double");
