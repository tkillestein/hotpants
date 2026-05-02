===========
Quick Start
===========

Basic Usage
===========

The simplest HOTPANTS command:

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim difference.fits

This produces a difference image where:

- ``science.fits`` is your science (observation) image
- ``template.fits`` is your reference (template) image
- ``difference.fits`` is the output difference image: D = I − T⊗K(x,y)

FFT-Accelerated Convolution
=============================

For large images, enable FFT acceleration (~3–8× speedup):

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim difference.fits -fft

Multi-Threading
================

Use multiple CPU cores (recommended):

.. code-block:: bash

   OMP_NUM_THREADS=4 hotpants -inim science.fits -tmplim template.fits -outim difference.fits -fft

Typical scaling: 2–4× faster on 4 cores. Auto-detects available cores if ``OMP_NUM_THREADS`` is unset.

Additional Outputs
===================

Save the noise map and convolved image:

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim difference.fits \
     -oni noise.fits -oci convolved.fits -fft

Output layers are stacked in a multi-layer FITS file. Common options:

- ``-oni FILE`` — Output noise map
- ``-oci FILE`` — Output convolved (template or science) image
- ``-ond FILE`` — Output noise-scaled difference
- ``-omi FILE`` — Output bad-pixel mask
- ``-allm`` — Include all available layers

Saving Kernel Information
==========================

Save the fitted kernel positions and values:

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim difference.fits \
     -savexy kernel_positions.txt -ki kernel_header.fits -fft

- ``-savexy FILE`` — Positions of stamps used for kernel fitting
- ``-ki FILE`` — Kernel coefficients in FITS header (can be reused with ``-ki`` on another image)

Processing Image Regions
=========================

For large images, divide into independent regions (better parallelization):

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim difference.fits \
     -nrx 2 -nry 2 -fft

This divides the image into a 2×2 grid of regions, each processed independently.

- ``-nrx N`` — Number of regions in X
- ``-nry N`` — Number of regions in Y

Masking Bad Pixels
===================

Use bad-pixel masks to exclude regions:

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim difference.fits \
     -imi science_mask.fits -tmi template_mask.fits -omi output_mask.fits -fft

- ``-imi FILE`` — Input image mask (0 = bad, non-zero = good)
- ``-tmi FILE`` — Template mask
- ``-omi FILE`` — Output bad-pixel mask

Noise Propagation
==================

If you have noise estimates (variance maps):

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim difference.fits \
     -ini science_noise.fits -tni template_noise.fits -oni output_noise.fits -fft

- ``-ini FILE`` — Input science noise map
- ``-tni FILE`` — Template noise map
- ``-oni FILE`` — Output noise map (propagated through convolution)

The output noise is computed as: σ²_out = σ²_I + (T⊗K)² ⊗ σ²_T

Common Scenarios
================

**Scenario 1: Basic image difference (default convolution)**

.. code-block:: bash

   hotpants -inim image.fits -tmplim template.fits -outim diff.fits

**Scenario 2: Large image, fast multi-threaded FFT**

.. code-block:: bash

   OMP_NUM_THREADS=8 hotpants -inim image.fits -tmplim template.fits -outim diff.fits -fft

**Scenario 3: Preserve science image PSF (sharp science, dull template)**

.. code-block:: bash

   hotpants -inim image.fits -tmplim template.fits -outim diff.fits -c i -fft

(Convolution is applied to template instead of science, preserving science PSF.)

**Scenario 4: Extract transients with noise scaling**

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim diff.fits \
     -ond diff_noise_scaled.fits -ni -ndm -fft

Output includes difference image, noise map, and noise-scaled difference.

**Scenario 5: Multiple regions with diagnostics**

.. code-block:: bash

   hotpants -inim image.fits -tmplim template.fits -outim diff.fits \
     -nrx 3 -nry 3 -savexy kernel_pos.txt -fft

Divides image into 3×3 regions, saves kernel fitting positions.

Next Steps
==========

- See :doc:`psf_matching` for guidance on when to convolve which image
- See :doc:`tuning` for performance optimization and parameter tuning
- See :doc:`../api/c_api` for full C API reference
