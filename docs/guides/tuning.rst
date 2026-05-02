=======================
Tuning & Performance
=======================

This guide covers parameter tuning for accuracy and performance optimization.

Kernel Fitting Parameters
==========================

``-ko`` (Kernel spatial order)
-------------------------------

Polynomial degree for spatial variation of kernel coefficients.

- Default: ``-ko 2`` (quadratic)
- Range: 0–6
- Effect: Higher = more spatial detail, but slower and prone to overfitting

Recommendation:

- For small images (< 2048²): ``-ko 2`` (default)
- For large images with smooth PSF: ``-ko 1`` (linear, faster)
- For complex PSF variation: ``-ko 3`` (cubic, more parameters)

Example:

.. code-block:: bash

   hotpants -inim image.fits -tmplim template.fits -outim diff.fits -ko 1 -fft  # Fast

``-bgo`` (Background spatial order)
-------------------------------------

Polynomial degree for background (DC level) variation within regions.

- Default: ``-bgo 1`` (linear)
- Range: 0–6
- Effect: Similar to ``-ko``; usually not the bottleneck

Recommendation: Keep at default (``-bgo 1``).

``-r`` (Kernel half-width)
----------------------------

Size of the convolution kernel in pixels.

- Default: ``-r 10`` (21×21 kernel)
- Range: 2–50 (larger kernels slow down fitting)
- Effect: Larger kernels capture broader PSF variations

Recommendation:

- Measure your PSF FWHM
- Set ``r ≈ 3–4 × FWHM``
- For FWHM = 3 px: ``-r 9``
- For FWHM = 6 px: ``-r 20``

.. code-block:: bash

   hotpants ... -r 15 -fft

``-ng`` (Kernel basis Gaussians)
---------------------------------

Number of Gaussian basis functions and their widths.

- Default: ``-ng 3 6 0.7 4 1.5 2 3.0`` (3 Gaussians: σ=0.7, 1.5, 3.0 px)
- Format: ``-ng count poly_deg_1 sigma_1 poly_deg_2 sigma_2 ...``

Recommendation:

- For broad kernels: Use your ``sigma_match`` as the middle Gaussian (see :doc:`psf_matching`)
- For narrow kernels: Use finer resolution (more basis functions)

Example (broad kernel):

.. code-block:: bash

   sigma_match=2.5
   hotpants ... -ng 3 6 0.5*$sigma_match 4 $sigma_match 2 2.0*$sigma_match -fft

Example (narrow kernel):

.. code-block:: bash

   hotpants ... -ng 4 6 0.5 5 1.0 4 1.5 3 2.5 -fft

Stamp Fitting & Rejection
==========================

``-ft`` (Fitting threshold)
-----------------------------

RMS threshold for centroid fitting in stamps.

- Default: ``-ft 20.0`` (sigma)
- Higher = fewer stamps used, cleaner fit, slower convergence
- Lower = more stamps used, noisier fit, risk of overfitting

Recommendation:

- For high-SNR data (clean images): ``-ft 15.0`` (aggressive)
- For noisy data (crowded fields): ``-ft 25.0`` (conservative)

.. code-block:: bash

   hotpants ... -ft 22.0 -fft

``-ks`` (Bad-stamp rejection threshold)
----------------------------------------

Sigma threshold for rejecting bad stamps during iterative fitting.

- Default: ``-ks 2.0`` (2-sigma clipping)
- Higher = more stamps kept, noisier fit
- Lower = fewer stamps kept, cleaner fit

Recommendation:

- Keep at default (``-ks 2.0``) for most cases
- Use ``-ks 3.0`` if you trust your data quality

``-nss`` (Number of substamps per stamp)
------------------------------------------

Number of bright stars to use per stamp for fitting.

- Default: ``-nss 3``
- Range: 1–20
- Effect: More substamps = better averaging, but slower if stars are scarce

Recommendation:

- For crowded fields: ``-nss 5`` (more stars available)
- For sparse fields: ``-nss 2`` (fewer bright stars)

.. code-block:: bash

   hotpants ... -nss 4 -fft

Image Regions
==============

``-nrx``, ``-nry`` (Region grid)
---------------------------------

Divide the image into independent regions for parallel processing.

- Default: ``-nrx 1 -nry 1`` (single region)
- Range: 1–10 per dimension
- Effect: More regions = better parallelization, smaller independent fits

Recommendation:

- For large images (> 4096²): ``-nrx 2 -nry 2`` or ``3 -nry 3``
- For multi-CPU systems: Set regions so each CPU handles ~1 region

Examples:

.. code-block:: bash

   # Quad-core CPU, large image
   hotpants ... -nrx 2 -nry 2 -fft

   # 8-core CPU, very large image
   hotpants ... -nrx 2 -nry 2 -fft -c i  # 4 independent regions, 8 threads = 2 threads/region

Stamp Grid
===========

``-nsx``, ``-nsy`` (Stamps per region)
----------------------------------------

Grid of stamps within each region (used for kernel fitting).

- Default: ``-nsx 10 -nsy 10`` (100 stamps per region)
- Range: 5–20
- Effect: More stamps = better fit, slower fitting

Recommendation:

- Keep at default (10×10) for most cases
- Increase to 12×12 for complex PSF variation
- Decrease to 8×8 for speed if fitting time is bottleneck

.. code-block:: bash

   hotpants ... -nsx 8 -nsy 8 -fft  # Faster, less detail

Performance Tuning
====================

FFT Acceleration
-----------------

**Always use FFT for large images** (> 1024×1024):

.. code-block:: bash

   hotpants ... -fft

Expected speedups:

- 1024×1024: **3–4×** faster
- 4096×4096: **5–8×** faster
- 8192×8192: **8–12×** faster

Multi-Threading
-----------------

**Use all available CPU cores:**

.. code-block:: bash

   OMP_NUM_THREADS=8 hotpants ... -fft

Typical scaling (on 4-core CPU):

- 1 thread: baseline (100%)
- 2 threads: ~1.8× faster
- 4 threads: ~3.5× faster (expected 4×, but some overhead)

Profiling Tips
==============

**Time individual steps:**

.. code-block:: bash

   time OMP_NUM_THREADS=4 hotpants ... -fft

**Profile with Linux perf (advanced):**

.. code-block:: bash

   perf stat -e cycles,instructions,cache-misses \
     OMP_NUM_THREADS=4 hotpants ... -fft

Look for:

- IPC (instructions per cycle) > 1.5 = good
- Cache-miss rate < 5% = good

Memory Usage
============

HOTPANTS memory footprint scales with image size and region count.

- Per-thread overhead: ~100 MB
- Per-image: ~8 × (image_size in pixels) + overhead
- FFT path adds ~3–4× for temporary FFT buffers

Recommendation:

- For 4096×4096 image on 4-thread system: ~1–2 GB RAM needed
- If memory-limited: reduce ``-nrx`` / ``-nry`` (process fewer regions in parallel)

.. code-block:: bash

   OMP_NUM_THREADS=2 hotpants ...  # Use 2 threads instead of 4 to reduce RAM

Typical Tuning Workflows
=========================

**Scenario 1: Large image, speed-critical (e.g., real-time survey)**

.. code-block:: bash

   OMP_NUM_THREADS=8 hotpants -inim image.fits -tmplim template.fits \
     -outim diff.fits -nrx 2 -nry 2 -nsx 8 -nsy 8 -ko 1 -r 12 -fft

Trade-off: Speed over accuracy. Simpler kernel, fewer stamps.

**Scenario 2: Archival data, accuracy-critical**

.. code-block:: bash

   hotpants -inim image.fits -tmplim template.fits -outim diff.fits \
     -nrx 1 -nry 1 -nsx 10 -nsy 10 -ko 2 -r 15 -ng 3 6 0.5 4 1.0 2 2.0 -ft 25.0

Trade-off: Accuracy over speed. Full kernel, good coverage.

**Scenario 3: Interactive/iterative (finding right PSF match)**

.. code-block:: bash

   hotpants -inim image.fits -tmplim template.fits -outim diff.fits \
     -nrx 1 -nry 1 -nsx 6 -nsy 6 -ko 1 -r 10 -fft

Trade-off: Balance. Fewer stamps for fast feedback.

Benchmarking
============

Time a typical run and track results:

.. code-block:: bash

   # Full timing
   /usr/bin/time -v hotpants ... -fft 2>&1 | grep -E "User|System|Elapsed|Maximum"

   # Simple elapsed time
   time hotpants ... -fft

Compare before/after optimization changes.

Next Steps
==========

- See :doc:`psf_matching` for PSF-matching decisions
- See :doc:`../api/c_api` for full option reference
- File an issue with profiling results if performance is unexpectedly poor
