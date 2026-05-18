/*
 * hotpants_config.h — Thread-safe configuration struct
 *
 * Defines hotpants_config_t: a struct that bundles all configuration for a
 * single fit_kernel / spatial_convolve call.  Including this header in both
 * hotpants_wrapper.c and hotpants_api.h avoids duplicate stamp_struct
 * definitions while still exposing the type to Python ctypes.
 */

#ifndef HOTPANTS_CONFIG_H
#define HOTPANTS_CONFIG_H

/**
 * @brief Configuration struct for thread-safe kernel-fitting calls.
 *
 * Pass a pointer to the _ctx wrapper functions in hotpants_api.h.
 * Each _ctx function applies these fields to the internal C globals atomically
 * (under a mutex), calls the core algorithm, then reads computed outputs back.
 *
 * Input fields must be filled by the caller before passing to initKernelGlobals_ctx().
 * Derived fields (nCompKer, nComp, ...) are written by initKernelGlobals_ctx().
 */
typedef struct {
    /* --- Kernel geometry --- */
    int    hwKernel;          /**< Half-width of kernel region (px) */
    int    kerOrder;          /**< Spatial polynomial order for kernel coefficients */
    int    bgOrder;           /**< Spatial polynomial order for background */

    /* --- Basis type --- */
    int    iBasisType;        /**< 0 = Gaussian, 1 = delta function */
    double rDeltaRegularization; /**< Delta basis Laplacian penalty weight (>= 0) */
    int    useTPS;            /**< 1 = use TPS spatial variation, 0 = polynomial */
    double tpsSmoothing;      /**< TPS regularisation parameter */

    /* --- Stamp layout --- */
    int    nKSStamps;         /**< Number of kernel test substamps per stamp */
    int    hwKSStamp;         /**< Half-width of kernel test substamp (px) */
    int    nRegX;             /**< Number of regions along x axis */
    int    nRegY;             /**< Number of regions along y axis */
    int    nStampX;           /**< Stamps per region along x */
    int    nStampY;           /**< Stamps per region along y */
    int    useFullSS;         /**< 1 = auto-divide region into kernel-sized substamps */

    /* --- Image thresholds --- */
    float  tUThresh;          /**< Template upper saturation threshold */
    float  tLThresh;          /**< Template lower valid data threshold */
    float  iUThresh;          /**< Science image upper saturation threshold */
    float  iLThresh;          /**< Science image lower valid data threshold */

    /* --- Noise model --- */
    float  tGain;             /**< Template gain (e-/ADU) */
    float  iGain;             /**< Science image gain (e-/ADU) */
    float  tRdnoise;          /**< Template read noise (e-) */
    float  iRdnoise;          /**< Science image read noise (e-) */
    float  tPedestal;         /**< Template baseline offset (ADU) */
    float  iPedestal;         /**< Science image baseline offset (ADU) */

    /* --- Fit quality --- */
    float  kerFitThresh;      /**< Sigma threshold for stamps used in fit */
    float  scaleFitThresh;    /**< Scale factor applied to fit threshold on retry */

    /* --- Execution control --- */
    int    verbose;           /**< Verbosity level (0 = silent) */
    int    nThread;           /**< Number of OpenMP threads */

    /* --- Derived fields (written by initKernelGlobals_ctx) --- */
    int    nCompKer;          /**< Number of kernel basis components */
    int    nComp;             /**< Number of spatial polynomial terms */
    int    nCompBG;           /**< Number of background polynomial terms */
    int    nCompTotal;        /**< Total component count */
    int    fwKernel;          /**< Full kernel width = 2*hwKernel + 1 */
    int    fwKSStamp;         /**< Full substamp width = 2*hwKSStamp + 1 */
    int    fwStamp;           /**< Stamp frame width (px) */
    int    kcStep;            /**< Convolution step size */

    /* --- Valid stamp count (set by build phase, consumed by fit phase) --- */
    int    nS;                /**< Number of valid stamps after buildStamps */
} hotpants_config_t;

#endif /* HOTPANTS_CONFIG_H */
