/**
 * @file basis.h
 * @brief Kernel basis function abstraction for pluggable kernel parameterizations.
 *
 * Provides a unified interface for different kernel representations:
 * - Gaussian basis (multi-Gaussian with polynomial spatial variation)
 * - Delta basis (pixel-level delta functions, Bramich 2008)
 * - Future: PCA basis, Fourier basis, etc.
 *
 * Each basis implements a set of operations:
 * 1. convolve_stamp() - correlate image with a basis function
 * 2. eval_kernel() - assemble kernel at image position from coefficients
 * 3. init() - initialize basis-specific state
 * 4. cleanup() - free basis-specific state
 */

#ifndef HOTPANTS_BASIS_H
#define HOTPANTS_BASIS_H

/* =====================================================================
   KERNEL BASIS ABSTRACTION
   ===================================================================== */

/**
 * @brief Kernel basis function interface.
 *
 * Each basis type (Gaussian, Delta, PCA, etc.) implements these operations
 * to enable generic matrix building and kernel evaluation.
 */
typedef struct kernel_basis {
  const char* name;           /**< Basis name: "gaussian", "delta", etc. */
  int nbasis;                 /**< Number of basis functions (nCompKer) */

  /**
   * @brief Convolve a stamp with a single basis function.
   *
   * @param[in] stamp Stamp pixel array (mStampX × mStampY)
   * @param[in] mStampX, mStampY Stamp dimensions
   * @param[in] basisIdx Index of basis function (0 to nbasis-1)
   * @param[out] pixFitted Output correlation array (mStampX × mStampY)
   * @return 0 on success, -1 on error
   *
   * For Gaussian basis: separable convolution with precomputed kernel_vec[]
   * For Delta basis: direct pixel extraction
   */
  int (*convolve_stamp)(double* stamp, int mStampX, int mStampY,
                        int basisIdx, double* pixFitted);

  /**
   * @brief Evaluate the spatial kernel at image position (xi, yi).
   *
   * @param[in] xi, yi Image coordinates
   * @param[in] kernelSol Fitted kernel coefficients (output of fitKernel)
   * @return Kernel pixel sum (flux-scaling factor)
   *
   * Assembles kernel[] global array from kernelSol and returns its sum.
   * Handles spatial variation (polynomial or TPS) internally.
   */
  double (*eval_kernel)(int xi, int yi, double* kernelSol);

  /**
   * @brief Initialize basis-specific state.
   *
   * Called after hwKernel is set but before fitKernel().
   * @return Number of basis functions (nbasis) on success, -1 on error
   */
  int (*init)(void);

  /**
   * @brief Clean up basis-specific resources.
   *
   * Called when switching bases or at program exit.
   * @return 0 on success
   */
  int (*cleanup)(void);

} kernel_basis_t;

/* =====================================================================
   BASIS REGISTRY AND DISPATCH
   ===================================================================== */

/** Global active basis (set by fitKernel based on iBasisType) */
extern kernel_basis_t* active_basis;

/**
 * @brief Get basis implementation for given type.
 *
 * @param[in] basis_type BASIS_TYPE_GAUSSIAN (0) or BASIS_TYPE_DELTA (1)
 * @return Pointer to basis implementation, or NULL if invalid
 */
kernel_basis_t* get_basis_for_type(int basis_type);

/**
 * @brief Set the active basis and initialize it.
 *
 * Calls old basis cleanup() if active, then calls new basis init().
 *
 * @param[in] basis_type BASIS_TYPE_GAUSSIAN (0) or BASIS_TYPE_DELTA (1)
 * @return 0 on success, -1 on error
 */
int set_active_basis(int basis_type);

/**
 * @brief Clean up the active basis.
 *
 * @return 0 on success
 */
int cleanup_active_basis(void);

/* =====================================================================
   BASIS IMPLEMENTATIONS (defined in basis_gaussian.c, basis_delta.c)
   ===================================================================== */

extern kernel_basis_t gaussian_basis;
extern kernel_basis_t delta_basis;

#endif  /* HOTPANTS_BASIS_H */
