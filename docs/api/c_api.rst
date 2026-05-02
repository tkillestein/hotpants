====
C API
====

Core Library Reference
=======================

The HOTPANTS core library provides functions for kernel fitting and spatially-varying convolution.

.. note::

   This documentation is auto-generated from Doxygen comments in the C source code using Breathe.

Main Functions
==============

Kernel Fitting
---------------

.. doxygenfunction:: fitKernel
   :project: HOTPANTS

Convolution (Direct)
---------------------

.. doxygenfunction:: spatial_convolve
   :project: HOTPANTS

Convolution (FFT-accelerated)
------------------------------

.. doxygenfunction:: spatial_convolve_fft
   :project: HOTPANTS

Helper Functions
-----------------

.. doxygenfunction:: make_kernel
   :project: HOTPANTS

.. doxygenfunction:: xy_conv_stamp
   :project: HOTPANTS

.. doxygenfunction:: getKernelVec
   :project: HOTPANTS

PSF Analysis
============

.. doxygenfunction:: getPsfCenters
   :project: HOTPANTS

.. doxygenfunction:: buildStamps
   :project: HOTPANTS

Statistics & Masking
====================

.. doxygenfunction:: getStampStats3
   :project: HOTPANTS

.. doxygenfunction:: makeNoiseImage4
   :project: HOTPANTS

Linear Algebra (LAPACK-based)
==============================

.. doxygenfunction:: build_matrix
   :project: HOTPANTS

.. doxygenfunction:: build_scprod
   :project: HOTPANTS

Legacy Functions (Numerical Recipes LU — deprecated)
======================================================

These functions are retained for compatibility but have been replaced by LAPACK Cholesky:

.. doxygenfunction:: ludcmp
   :project: HOTPANTS

.. doxygenfunction:: lubksb
   :project: HOTPANTS

All Files
=========

.. toctree::
   :maxdepth: 1

.. doxygenindex::
   :project: HOTPANTS

See Also
========

- :doc:`../guides/algorithm` — Mathematical overview
- :doc:`../guides/tuning` — Parameter reference for kernel fitting
- `GitHub Repository <https://github.com/tkillestein/hotpants>`_
