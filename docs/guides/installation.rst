============
Installation
============

HOTPANTS is built with CMake.

CMake
=====

**Prerequisites:**

- CMake ≥ 3.10
- C compiler with C99 support (GCC, Clang)
- CFITSIO ≥ 3.45 (FITS I/O)
- OpenBLAS or LAPACK (linear algebra)
- FFTW3 ≥ 3.3 (optional, for FFT acceleration)
- OpenMP (optional, for multi-threading)

**Linux (Ubuntu/Debian):**

.. code-block:: bash

   sudo apt-get install cmake build-essential libcfitsio-dev libopenblas-dev libfftw3-dev libomp-dev

**macOS (Homebrew):**

.. code-block:: bash

   brew install cmake cfitsio openblas fftw open-mpi

**Build:**

.. code-block:: bash

   cmake -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build
   ./build/hotpants -h  # Verify installation

**Optional: Install to system:**

.. code-block:: bash

   sudo cmake --install build
   hotpants -h  # Should be in PATH

Configuration options:

- ``-DUSE_FFTW=ON/OFF`` — Enable FFT acceleration (default: ON if FFTW3 found)
- ``-DUSE_OPENMP=ON/OFF`` — Enable multi-threading (default: ON if found)
- ``-DCMAKE_C_FLAGS="-O3 -march=native -funroll-loops"`` — Compiler optimization

Docker (Future)
===============

A Dockerfile is planned for containerized builds. Check back later or file an issue to request it.

Verifying Installation
======================

After installation, verify HOTPANTS works:

.. code-block:: bash

   hotpants -h

Should display the full option reference.

Quick sanity check:

.. code-block:: bash

   # Create test FITS images (requires astropy)
   python3 << 'EOF'
   import numpy as np
   from astropy.io import fits

   # Create simple test images
   data = np.random.normal(1000, 10, (512, 512)).astype(np.float32)

   fits.writeto('science.fits', data, overwrite=True)
   fits.writeto('template.fits', data + 5, overwrite=True)
   EOF

   # Run HOTPANTS
   hotpants -inim science.fits -tmplim template.fits -outim diff.fits -fft

   # Check output
   python3 -c "from astropy.io import fits; print(fits.getdata('diff.fits').shape)"

Troubleshooting
===============

**"command not found: hotpants"**

Make sure the build directory is in ``PATH``, or use the full path: ``./build/hotpants``.

**"CFITSIO not found"**

Install CFITSIO and ensure CMake can find it:

.. code-block:: bash

   pkg-config --cflags --libs cfitsio  # Should output paths

If not found, edit ``CMakeLists.txt`` to specify explicit paths.

**"Doxygen not found" (when building docs)**

Install Doxygen: ``brew install doxygen`` (macOS) or ``sudo apt-get install doxygen`` (Linux).

Next Steps
==========

- See :doc:`quickstart` for usage examples
- See :doc:`psf_matching` for PSF matching strategies
- See :doc:`tuning` for performance optimization
