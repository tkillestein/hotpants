typedef struct
{
   int       x0,y0;       /* origin of stamp in region coords*/
   int       x,y;         /* center of stamp in region coords*/
   int       nx,ny;       /* size of stamp */
   int       *xss;        /* x location of test substamp centers */
   int       *yss;        /* y location of test substamp centers */
   int       nss;         /* number of detected substamps, 1 .. nss     */
   int       sscnt;       /* represents which nss to use,  0 .. nss - 1 */
   double    **vectors;   /* contains convolved image data */
   double    *krefArea;   /* contains kernel substamp data */
   double    **mat;       /* fitting matrices */
   double    *scprod;     /* kernel sum solution */
   double    sum;         /* sum of fabs, for sigma use in check_stamps */
   double    mean;
   double    median;
   double    mode;        /* sky estimate */
   double    sd;
   double    fwhm;
   double    lfwhm;
   double    chi2;        /* residual in kernel fitting */
   double    norm;        /* kernel sum */
   double    diff;        /* (norm - mean_ksum) * sqrt(sum) */
} stamp_struct;

/* Include main.c defines HOTPANTS_DEFINE_GLOBALS before including this header
   so that variables are defined once there and declared extern everywhere else. */
#ifdef HOTPANTS_DEFINE_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

/* GLOBAL VARS POSSIBLY SET ON COMMAND LINE */
EXTERN char      *template, *image, *outim;

EXTERN float     tUThresh, tUKThresh, tLThresh, tGain, tRdnoise, iUThresh, iUKThresh, iLThresh, iGain, iRdnoise;
EXTERN char      *tNoiseIm, *iNoiseIm, *tMaskIm, *iMaskIm, *kernelImIn, *kernelImOut, *outMask;
EXTERN float     tPedestal, iPedestal;
EXTERN int       hwKernel;
EXTERN float     kerFitThresh, scaleFitThresh, minFracGoodStamps;
EXTERN float     kfSpreadMask1, kfSpreadMask2;
EXTERN int       gdXmin, gdXmax, gdYmin, gdYmax;
EXTERN int       nRegX, nRegY;
EXTERN char      *regFile;
EXTERN char      *regKeyWord;
EXTERN int       numRegKeyWord;
EXTERN int       nStampY, nStampX, useFullSS;
EXTERN int       nKSStamps, hwKSStamp;
EXTERN char      *sstampFile;
EXTERN int       findSSC;
EXTERN int       kerOrder, bgOrder;
EXTERN float     statSig, kerSigReject, kerFracMask;
EXTERN char      *forceConvolve, *photNormalize, *figMerit;
EXTERN int       sameConv, rescaleOK;
EXTERN float     fillVal, fillValNoise;
EXTERN char      *effFile, *noiseImage, *sigmaImage, *convImage;
EXTERN int       doSum, inclNoiseImage, inclSigmaImage, inclConvImage, noClobber;
EXTERN int       doKerInfo, outShort, outNShort;
EXTERN float     outBzero, outBscale, outNiBzero, outNiBscale;
EXTERN int       convolveVariance;
EXTERN int       usePCA, fwKernelPCA;
EXTERN float     **PCA;

/* GLOBAL VARS NOT SET ON COMMAND LINE */
EXTERN int       ngauss, *deg_fixe;
EXTERN float     *sigma_gauss;

EXTERN int       rPixX, rPixY;
EXTERN int       nStamps, nS, nCompKer, nC;

EXTERN int       nComp, nCompBG, nBGVectors, nCompTotal;

EXTERN int       fwKernel, fwStamp, hwStamp, fwKSStamp, kcStep;
EXTERN int       cmpFile;
EXTERN float     *temp, *temp2;
EXTERN double    *check_stack,*filter_x,*filter_y,**kernel_vec;
EXTERN double    **wxy,*kernel_coeffs,*kernel,**check_mat,*check_vec;
EXTERN char      version[32];

/* REGION SIZED */
EXTERN int       *mRData;   /* bad input data mask */

/* armin */
/* a dummy varialbe to do some testing */
EXTERN int        dummy;
/* verbose for debugging */
EXTERN int        verbose;
/* cmp file stuff */
EXTERN char       xyfilename[1000];
EXTERN int        savexyflag;
EXTERN float      *xcmp,*ycmp;
EXTERN int        Ncmp;
