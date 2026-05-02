============================================
HOTPANTS — High Order Transform of PSF ANd Template Subtraction
============================================

A fast, production-grade image differencing tool for astronomical image subtraction.

Welcome to the HOTPANTS documentation. Whether you're getting started or tuning advanced parameters, you'll find comprehensive guides and API references here.

.. toctree::
   :maxdepth: 2
   :caption: User Guides

   guides/installation
   guides/quickstart
   guides/psf_matching
   guides/tuning

.. toctree::
   :maxdepth: 2
   :caption: Reference

   guides/algorithm
   api/c_api
   api/python_api

.. toctree::
   :maxdepth: 1
   :caption: Contributing

   Contributing (future) <https://github.com/tkillestein/hotpants/blob/main/CONTRIBUTING.md>
   Development Guide <https://github.com/tkillestein/hotpants/blob/main/CLAUDE.md>

---

Quick Links
===========

- **GitHub Repository**: https://github.com/tkillestein/hotpants
- **Issue Tracker**: https://github.com/tkillestein/hotpants/issues
- **ASCL Entry**: https://ascl.net/1504.004
- **Citation**: `Alard & Lupton (1998) <https://iopscience.iop.org/article/10.1086/305984>`_

Installation
============

Install via CMake:

.. code-block:: bash

   cmake -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build
   ./build/hotpants [options]

See :doc:`guides/installation` for detailed instructions and dependencies.

Quick Start
===========

.. code-block:: bash

   # Basic usage
   hotpants -inim science.fits -tmplim template.fits -outim difference.fits

   # With FFT acceleration (much faster)
   hotpants -inim science.fits -tmplim template.fits -outim difference.fits -fft

   # Multi-threaded
   OMP_NUM_THREADS=4 hotpants -inim science.fits -tmplim template.fits -outim difference.fits -fft

See :doc:`guides/quickstart` for more examples.

---

API Reference
=============

.. note::
   C API documentation is auto-generated from source code comments using Doxygen and Breathe.

The C core library provides:

- **Kernel fitting** (``fitKernel()``) — Least-squares polynomial kernel fit per region
- **Spatial convolution** (``spatial_convolve()``, ``spatial_convolve_fft()``) — Apply fitted kernel to images
- **PSF analysis** (``getPsfCenters()``) — Locate bright stars for kernel fitting
- **Image I/O** — FITS read/write via CFITSIO

See :doc:`api/c_api` for full reference.

Python bindings (coming soon) will wrap the C core for numpy array inputs.

---

Algorithm
=========

HOTPANTS solves:

.. math::

   D = I - T \otimes K(x, y)

where:

- D is the difference image
- I is the science image
- T is the template image
- K(x, y) is a spatially-varying convolution kernel

The kernel is parameterized as a sum of Gaussian basis functions with spatial polynomial coefficients. See :doc:`guides/algorithm` for details.

---

Performance
===========

On a modern 4-core CPU:

- **1024×1024 image**: ~0.1 s (FFT-accelerated)
- **4096×4096 image**: ~2 s (FFT-accelerated)

FFT acceleration provides **3–8× speedup** over direct convolution. See :doc:`guides/tuning` for performance optimization.

---

License
=======

MIT License (2013–2026). See LICENSE in the repository.
