==========
Python API
==========

Python Bindings (Coming Soon)
==============================

HOTPANTS will provide a modern Python API wrapping the C core library for direct numpy array processing.

**Status**: Under development. See `CLAUDE.md <https://github.com/tkillestein/hotpants/blob/main/CLAUDE.md>`_ for roadmap.

Planned API
===========

Once available, the Python interface will look like:

.. code-block:: python

   import hotpants
   import numpy as np

   # Load images
   science = fits.getdata('science.fits')
   template = fits.getdata('template.fits')

   # Configure parameters
   params = hotpants.HOTPANTSParams(
       kernel_half_width=10,
       regions_x=2,
       regions_y=2,
       use_fft=True,
   )

   # Run subtraction
   diff, noise, convolved = hotpants.subtract(
       science, template, params
   )

   # Save results
   fits.writeto('diff.fits', diff)
   fits.writeto('noise.fits', noise)

Key Features:

- **Direct numpy array input/output** — no FITS I/O in the binding
- **Parameter dataclass** — clean configuration API
- **Backward compatibility** — CLI continues to work alongside Python API
- **Minimal overhead** — thin cffi or pybind11 layer

Implementation Details
======================

**Binding strategy**: `cffi <https://cffi.readthedocs.io/>`_ (preferred) or `pybind11 <https://pybind11.readthedocs.io/>`_

- **cffi**: Avoids Python/C ABI fragility; easier to maintain across Python versions
- **pybind11**: Requires C++ shim; cleaner C++ integration if needed

**Under the hood**:

1. User calls ``hotpants.subtract(I, T, params)`` with numpy arrays
2. Binding allocates FITS-like in-memory buffers
3. Calls C core ``main()`` logic (region loop, fitting, convolution)
4. Unmarshals results back to numpy arrays
5. Returns (difference_image, noise_map, ...)

See `pyproject.toml <https://github.com/tkillestein/hotpants/blob/main/pyproject.toml>`_ and `CLAUDE.md <https://github.com/tkillestein/hotpants/blob/main/CLAUDE.md>`_ for details.

Timeline
========

- **Current**: C CLI fully functional and optimized
- **Next phase**: Python bindings (cffi/pybind11)
- **Future**: Upstream to conda-forge, pip

Contributing
=============

Interested in helping implement Python bindings? See `CONTRIBUTING.md <https://github.com/tkillestein/hotpants/blob/main/CONTRIBUTING.md>`_ (coming soon) or file an issue on GitHub.

Related
=======

- :doc:`c_api` — C library reference
- :doc:`../guides/algorithm` — Algorithm details
- `Alard & Lupton (1998) <https://iopscience.iop.org/article/10.1086/305984>`_ — Original paper
