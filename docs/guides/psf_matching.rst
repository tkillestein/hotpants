=======================
PSF Matching Strategy
=======================

One of the most important choices in HOTPANTS is which image to convolve. This guide walks through the decision.

Core Problem
============

HOTPANTS solves: **D = I − T⊗K(x,y)**

where:

- **I** = science (observed) image
- **T** = template (reference) image
- **K(x,y)** = spatially-varying convolution kernel

The kernel K must match the science and template PSFs. Depending on their relative sizes, different strategies are needed.

Decision Tree
=============

**Step 1: Measure PSF widths**

Measure the PSF standard deviations (σ) in both images. Common tools:

- `astropy.photometry.DAOStarFinder` (Python)
- `SExtractor <https://www.astromatic.net/software/sextractor/>`_
- `PSFEx <https://www.astromatic.net/software/psfex/>`_

Example (Python):

.. code-block:: python

   from astropy.io import fits
   from photutils.detection import DAOStarFinder

   data = fits.getdata('science.fits')
   dao = DAOStarFinder(threshold=5*np.std(data))
   sources = dao(data)
   print(f"Median FWHM: {np.median(sources['fwhm'])}")  # Convert to σ if needed

**Step 2: Compare σ_science vs σ_template**

+-----------------------------+---------------------------------------------+
| Condition                   | Action                                      |
+=============================+=============================================+
| σ_science ≈ σ_template      | ✓ Proceed with default; any choice works    |
+-----------------------------+---------------------------------------------+
| σ_science < σ_template      | ⚠ Science sharper; see "Sharpened Science"  |
+-----------------------------+---------------------------------------------+
| σ_science > σ_template      | ✓ Ideal case; see "Natural Matching"        |
+-----------------------------+---------------------------------------------+

Case 1: Natural Matching (σ_science > σ_template)
===================================================

**Best case scenario.** The science image has a broader PSF than the template.

The kernel will *convolve the template* to match the science image's PSF. This is physically well-behaved: convolution is smoothing, which is stable.

**Recommended setup:**

1. Measure the matching Gaussian:

   .. math::

      \sigma_{\text{match}} = \sqrt{\sigma_I^2 - \sigma_T^2}

2. Set kernel basis Gaussians centered on σ_match:

   .. code-block:: bash

      hotpants -inim science.fits -tmplim template.fits -outim diff.fits \
        -ng 3 6 \
        0.5*$sigma_match 4 \
        $sigma_match 2 \
        2.0*$sigma_match -fft

   Where:

   - ``-ng 3 6`` = 3 basis functions, degree-6 polynomials
   - Arguments are: σ₁, degree₁, σ₂, degree₂, σ₃, degree₃
   - σ₁ = 0.5 × σ_match (finest detail)
   - σ₂ = σ_match (central Gaussian)
   - σ₃ = 2.0 × σ_match (broadest smoothing)

3. Inspect results:

   .. code-block:: bash

      fitsio -r diff.fits | head -20  # Check for artifacts

Case 2: Sharpened Science (σ_science < σ_template)
====================================================

**Problematic scenario.** The science image is sharper than the template, requiring **deconvolution** (sharpening).

Deconvolution is numerically unstable and prone to ringing artifacts and noise amplification.

**Option A: Convolve science instead (recommended)**

Force HOTPANTS to convolve the science image, avoiding deconvolution:

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim diff.fits -c i -fft

This solves: **D = I⊗K_c − T** where K_c is fitted to sharpen the template.

*Caveat:* This modifies the science image pixels. Use if the template is better-characterized.

**Option B: Pre-process science image**

Before matching, convolve the science image with its own PSF to make it broader:

.. code-block:: python

   from scipy.ndimage import convolve
   from astropy.io import fits

   science = fits.getdata('science.fits')

   # Create Gaussian kernel of science PSF
   sigma_science = 2.0  # pixels
   kernel = create_gaussian_kernel(sigma_science)

   # Convolve science with its own PSF
   science_convolved = convolve(science, kernel)

   # Now σ_new ≈ sqrt(2) × σ_science, often > σ_template
   fits.writeto('science_preconvolved.fits', science_convolved, overwrite=True)

Then run HOTPANTS normally:

.. code-block:: bash

   hotpants -inim science_preconvolved.fits -tmplim template.fits -outim diff.fits -fft

This approach is cleaner: the science image is modified *before* HOTPANTS, and the kernel fit is stable.

**Option C: Lower fitting threshold (careful)**

If deconvolution is unavoidable, use a higher fitting threshold to reject noisy stamps:

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim diff.fits \
     -ft 30.0 -ks 3.0 -fft

Parameters:

- ``-ft 30.0`` — RMS threshold for fitting (default 20.0); higher = fewer stamps used
- ``-ks 3.0`` — Sigma rejection for bad stamps (default 2.0); more aggressive clipping

**Use with caution:** Fewer stamps means worse spatial variation of the kernel.

Case 3: Identical PSFs (σ_science ≈ σ_template)
=================================================

In this ideal case, K ≈ δ (delta function) — no real convolution needed.

Run HOTPANTS normally:

.. code-block:: bash

   hotpants -inim science.fits -tmplim template.fits -outim diff.fits -fft

The difference image should be clean with minimal artifacts.

Kernel Basis Selection
=======================

For Cases 1 & 3, choose the kernel basis Gaussians:

.. code-block:: bash

   hotpants ... -ng 3 6 sigma_1 4 sigma_2 2 sigma_3

Default (for broad σ_match):

- σ₁ = 0.5 × σ_match (fine detail)
- σ₂ = 1.0 × σ_match (central, 6th-order polynomials)
- σ₃ = 2.0 × σ_match (broad, 2nd-order polynomials)

For sharper kernels (Case 2, Option A with ``-c i``), use finer basis:

.. code-block:: bash

   hotpants ... -c i -ng 3 6 0.3*$sigma_match 4 0.7*$sigma_match 2 1.5*$sigma_match

Diagnosing Problems
====================

**Symptom: Ringing or checkerboard patterns in difference image**

- Usually indicates overfitting or unstable kernel fit
- Reduce polynomial degree: ``-ko 1`` (instead of default 2)
- Increase fitting threshold: ``-ft 25.0``
- Check for cosmic rays in stamps: use ``-imi`` mask

**Symptom: Statistically correlated residuals**

- Kernel basis may be mistuned
- Try different σ values; center on actual σ_match
- Increase number of regions: ``-nrx 2 -nry 2``

**Symptom: Science image appears dimmer in difference**

- Indicates incomplete PSF matching; kernel is too weak
- Check that your PSF measurements are accurate
- Try manually increasing σ basis: ``2.0*sigma_match 4 3.0*sigma_match 2 4.0*sigma_match``

Validation Checklist
====================

After subtraction, verify quality:

1. **Visual inspection**: Does the difference image look smooth? Any obvious artifacts?
2. **Source statistics**: Are faint stars suppressed to noise level? Any residual halos?
3. **Noise scaling**: Is the noise level consistent with σ_I² + σ_T² (propagated)?
4. **Spatial variation**: Do kernel coefficients vary smoothly across the image?

.. code-block:: bash

   # Check noise in difference image
   python3 << 'EOF'
   from astropy.io import fits
   import numpy as np

   diff = fits.getdata('diff.fits')
   noise = fits.getdata('diff_noise.fits')

   # Measure RMS in clean region (away from sources)
   y, x = diff.shape
   clean = diff[y//4:3*y//4, x//4:3*x//4]
   rms_measured = np.std(clean[~np.isnan(clean)])
   rms_expected = np.sqrt(np.mean(noise**2))

   print(f"Measured RMS: {rms_measured:.3f}")
   print(f"Expected RMS: {rms_expected:.3f}")
   print(f"Ratio: {rms_measured/rms_expected:.2f}")
   EOF

Next Steps
==========

- See :doc:`tuning` for kernel parameter optimization
- See :doc:`../api/c_api` for detailed option reference
- File an issue on GitHub if your specific case is not covered
