==================
Algorithm Overview
==================

This document provides a mathematical overview of the HOTPANTS algorithm. For implementation details, see the C API documentation.

Problem Statement
=================

Given:

- **I** — science (observed) image
- **T** — template (reference) image

Find a kernel **K(x,y)** such that the difference image

.. math::

   D(x, y) = I(x, y) - [T \otimes K](x, y)

is minimized in the least-squares sense, where ⊗ denotes 2D convolution.

The kernel is spatially-varying: K changes across the image to account for PSF variations (e.g., from atmospheric turbulence or instrumental effects).

Kernel Parameterization
========================

Rather than fitting K at every pixel (infeasible), HOTPANTS models the kernel as a sum of **fixed Gaussian basis functions** with spatially-varying polynomial coefficients:

.. math::

   K(x, y) = \sum_{i=1}^{N_b} c_i(x, y) \cdot \phi_i

where:

- **φᵢ** — fixed Gaussian basis functions (e.g., σ = 0.7, 1.5, 3.0 px)
- **cᵢ(x,y)** — low-order polynomials in image coordinates (default: degree 6, 4, 2)
- **Nₑ** — number of basis functions (default: 3)

This reduces a 2D image-sized problem to fitting ~30–50 polynomial coefficients across the image.

Pipeline
=========

**Step 1: Divide into regions**

Split the image into rectangular regions (e.g., 2×2 grid for ``-nrx 2 -nry 2``) to handle spatially-varying PSFs locally.

**Step 2: Lay stamp grid**

Within each region, lay a grid of stamps (typically 10×10) where kernel fitting will occur.

**Step 3: Find bright stars**

In each stamp, identify bright point sources using centroid detection. These become the fitting points because:

- Point sources are relatively unaffected by the template PSF
- Their positions are determined precisely
- They provide good signal for constraining K

Function: ``getPsfCenters()`` in ``functions.c``.

**Step 4: Build normal equations**

For each stamp, convolve the template T with each basis element φᵢ:

.. math::

   T_i = T \otimes \phi_i

and stack the convolutions horizontally. For a stamp with bright stars at positions {xⱼ, yⱼ}, solve the least-squares problem:

.. math::

   c = \arg\min_c \sum_{\text{stamps}} \sum_{j} \left[ I(x_j, y_j) - \sum_i c_i(x_j, y_j) T_i(x_j, y_j) \right]^2

This accumulates into a symmetric positive-definite normal equation matrix:

.. math::

   \mathbf{A} \mathbf{c} = \mathbf{b}

where:

- **A** = ΦᵀΦ (design matrix)
- **b** = Φᵀ I (data vector)
- **Φ** = [φ₁, φ₂, ..., φₙᵦ] (basis matrix)

Functions: ``build_matrix()`` and ``build_scprod()`` in ``alard.c``.

**Step 5: Solve for coefficients**

Solve **Ac = b** using **Cholesky decomposition** (since A is positive-definite):

.. math::

   \mathbf{A} = \mathbf{L} \mathbf{L}^\top
   \mathbf{L} \mathbf{L}^\top \mathbf{c} = \mathbf{b}

Two triangular solves: Lx = b, then L^T c = x.

Function: ``fitKernel()`` in ``alard.c`` (calls LAPACK ``dpotrf`` + ``dpotrs``).

**Step 6: Sigma-clip bad stamps (iterate)**

Stamps with large residuals are rejected and fitting is repeated, typically 2–3 iterations.

**Step 7: Evaluate kernel across full image**

Using the fitted polynomial coefficients cᵢ, evaluate K(x,y) at every pixel:

.. math::

   K(x, y) = \sum_i c_i(x, y) \cdot \phi_i(x - x_0, y - y_0)

where φᵢ is re-centered at (x₀, y₀) for efficient computation.

**Step 8: Apply spatially-varying convolution**

Compute the difference image:

.. math::

   D(x, y) = I(x, y) - \sum_i c_i(x, y) [T \otimes \phi_i](x, y)

This is the slowest step (~60% CPU time in direct implementation).

Two implementations:

1. **Direct convolution** (``spatial_convolve()`` in ``alard.c``)
   - O(N k²) where N = image size, k = kernel half-width

2. **FFT acceleration** (``spatial_convolve_fft()`` in ``alard.c``)
   - Precomputes I⊗φᵢ for each basis via FFT (O(N log N))
   - Weighted sum: D = I − Σᵢ cᵢ(x,y) [I⊗φᵢ]
   - 3–8× faster for large images

**Step 9: Output**

Write difference image and optional layers (noise, convolved image, mask) to FITS.

Numerical Stability
====================

**PSF-matched kernels are well-conditioned** because:

- The normal equation matrix A = ΦᵀΦ is Gram matrix of basis functions
- Gaussian basis functions have controlled overlap and norm
- Condition number is typically κ(A) ≈ 10–100 (excellent)

**Cholesky decomposition is stable** for such matrices:

- Error amplification ≈ κ(A) × machine epsilon
- Single-precision (float) is usually sufficient; double precision used for robustness

**Deconvolution (when σ_I < σ_T) is ill-posed:**

- Requires sharpening the template (inverse of convolution)
- Noise amplified by factor ~(σ_T/δσ)²
- Not recommended; see :doc:`psf_matching` for alternatives

Noise Propagation
==================

If noise maps (variance) are provided via ``-ini`` and ``-tni``, the output noise is computed as:

.. math::

   \sigma_D^2 = \sigma_I^2 + (T \otimes K)^2 \otimes \sigma_T^2

where (T⊗K)² is element-wise squaring, and the convolution of variance follows from error propagation.

Implementation: ``makeNoiseImage4()`` in ``functions.c``.

Performance Bottlenecks
=======================

From profiling (see ``PROF`` in repo):

| Stage | CPU time | Notes |
|-------|----------|-------|
| ``spatial_convolve()`` direct | ~60% | Primary bottleneck; FFT helps here |
| ``getPsfCenters()`` | ~35% | Secondary; SIMD vectorization could help |
| ``xy_conv_stamp()`` | ~18% | Subset of PSF-centre finding |
| ``build_matrix()`` | ~8% | Fast; BLAS helps |
| Cholesky solve | < 5% | Negligible with LAPACK |

With FFT acceleration:

- ``spatial_convolve_fft()`` ~10–20% CPU (3–8× faster)
- ``getPsfCenters()`` remains ~35% (same code, no FFT benefit)

Future optimizations:

- SIMD vectorization of ``getPsfCenters()`` centroid loops (~50% of secondary bottleneck)
- Adaptive basis selection via PCA (reduce number of basis functions)
- GPU offload of stamp comparisons (modest benefit due to data transfer cost)

References
==========

.. [Alard1998] Alard, C. and Lupton, R. H., "A Method for Optimal Image Subtraction," *ApJ* **503**:325 (1998). `doi:10.1086/305984 <https://doi.org/10.1086/305984>`_

The original paper derives the least-squares kernel fitting and provides the foundational equations. HOTPANTS implements this with extensions for spatial variation (polynomial coefficients) and numerical optimization (Cholesky, FFT).

Next Steps
==========

- See :doc:`../api/c_api` for implementation details (function signatures, code-level invariants)
- See :doc:`tuning` for practical parameter choices
- See :doc:`psf_matching` for decision guidance on PSF matching
