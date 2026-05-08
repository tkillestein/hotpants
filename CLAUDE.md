# HOTPANTS — Development Guide

**High Order Transform of PSF ANd Template Subtraction**

HOTPANTS produces difference images between an astronomical science image and a
template by fitting and applying a spatially-varying convolution kernel that
matches the two PSFs. The core algorithm is from Alard & Lupton (1998), ApJ
503:325 (https://iopscience.iop.org/article/10.1086/305984).

License: MIT (Andy Becker, 2013). See `LICENSE`.

---

## Project Status (May 2026)

**Core implementation (C):**
- ✓ CMake build system with CFITSIO, OpenBLAS, FFTW3, OpenMP detection
- ✓ LAPACK Cholesky solver, OpenMP parallelism, FFT-accelerated convolution
- ✓ Doxygen-documented core functions, global naming legend, named constants
- ✓ Logging macros, memory allocation wrappers, contiguous allocation

**Python API (ctypes bindings):**
- ✓ C API header (`hotpants_api.h`) with function signatures and data structures
- ✓ Python module with `KernelConfig`, `RegionLayout`, `NoiseThresholds`, `KernelSolution` dataclasses
- ✓ Full ctypes binding wrapper (`_hotpants.py`) with numpy marshalling
- ✓ All unit tests passing (33/33), segfaults fixed, parameter validation working

**Documentation & build:**
- ✓ `pyproject.toml` configured for Python ≥ 3.11
- ✓ Sphinx + Breathe setup, Doxygen integration
- ✓ Test suite (`test_api.py`, `test_api_integration.py`, `test_regression.py`)

---

## Repository Structure

```
hotpants/
├── src/
│   ├── main.c              # Pipeline orchestration: FITS I/O, region loop, output assembly
│   ├── alard.c             # Core algorithm: kernel fitting, spatial convolution, LAPACK solve
│   ├── functions.c         # Stamps, PSF-centre finding, statistics, masking
│   ├── vargs.c             # CLI argument parsing (~50 options)
│   ├── maskim.c            # Standalone utility: apply mask to image
│   ├── globals.c           # Global variable definitions
│   ├── defaults.h          # Compile-time parameter defaults (with named constants)
│   ├── globals.h           # Global variable declarations + prefix legend
│   ├── functions.h         # Function prototypes and CFITSIO includes
│   └── hotpants_api.h      # Public C API for Python bindings (Doxygen-documented)
├── tests/
│   ├── conftest.py         # pytest configuration
│   ├── test_api.py         # Unit tests for Python API (KernelConfig, validation)
│   ├── test_api_integration.py  # End-to-end integration tests
│   ├── test_regression.py  # Regression tests against known outputs
│   └── __init__.py
├── docs/
│   ├── conf.py             # Sphinx configuration with Breathe (Doxygen integration)
│   ├── index.rst           # Documentation root
│   ├── api/                # API reference (generated from Doxygen + docstrings)
│   ├── guides/             # User guides and tutorials
│   ├── _static/            # Static assets (CSS, images)
│   └── _templates/         # Sphinx templates
├── CMakeLists.txt          # CMake build configuration
├── Doxyfile                # Doxygen configuration for C API documentation
├── pyproject.toml          # Python project metadata (uv-managed)
├── README.md               # User-facing option reference and quick start
├── CLAUDE.md               # This file: development guide for AI assistants
├── CONTRIBUTING.md         # (planned) Contributor guidelines
└── NOTES                   # Version history and changelog
```

### Key Functions

| Function                            | File          | Purpose                                                     | Status               |
|-------------------------------------|---------------|-------------------------------------------------------------|----------------------|
| `main()`                            | `main.c`      | Orchestrates the full pipeline                              | ✓                    |
| `fitKernel()`                       | `alard.c`     | Least-squares kernel fit per region                         | ✓ Documented         |
| `build_matrix()` / `build_scprod()` | `alard.c`     | Accumulate normal equations                                 | ✓ Documented         |
| `dpotrf`/`dpotrs` (LAPACK)          | `alard.c`     | Cholesky solve of normal equations                          | ✓ Replaced LU        |
| `spatial_convolve()`                | `alard.c`     | Direct convolution (primary bottleneck, ~60% CPU)           | ✓ Documented         |
| `spatial_convolve_fft()`            | `alard.c`     | FFT-accelerated convolution (3–8× faster)                   | ✓                    |
| `make_kernel()`                     | `alard.c`     | Evaluate kernel at image position from polynomials          | ✓ Documented         |
| `xy_conv_stamp()`                   | `alard.c`     | Separable 2D convolution with Gaussian basis                | ✓ Documented         |
| `getKernelVec()`                    | `alard.c`     | Initialise multi-Gaussian kernel basis vectors              | ✓ Documented         |
| `buildStamps()`                     | `functions.c` | Create stamp grid and compute statistics                    | ✓ In API             |
| `getPsfCenters()`                   | `functions.c` | Locate bright star centres (secondary bottleneck, ~35% CPU) | ✓ In API             |
| `getStampStats3()`                  | `functions.c` | Sigma-clipped statistics and histogram FWHM                 | ✓ In API, Documented |
| `makeNoiseImage4()`                 | `functions.c` | Propagate noise through convolution                         | ✓ Documented         |
| `vargs()`                           | `vargs.c`     | Parse all command-line options                              | ✓                    |

---

## Algorithm Overview

HOTPANTS solves **D = I − T⊗K(x,y)** where D is the difference image, I the
science image, T the template, and K(x,y) a spatially-varying convolution
kernel.

The kernel is parameterised as a sum of Gaussian basis functions weighted by
spatial polynomials:

```
K(x,y) = Σᵢ cᵢ(x,y) · φᵢ
```

where φᵢ are fixed Gaussian profiles (default: 3 Gaussians with σ = 0.7, 1.5,
3.0 pixels and polynomial degrees 6, 4, 2) and cᵢ(x,y) are low-order
polynomials in image position. This parameterization follows **Alard & Lupton
(1998), Section 2.2**.

**Pipeline steps:**

1. Divide image into regions (`-nrx`, `-nry`)
2. Within each region, lay a stamp grid (`-nsx`, `-nsy`)
3. In each stamp, locate bright PSF centres (`getPsfCenters`) — identifies ~3 substamps
   per stamp
4. For each substamp, convolve I (or T) with every basis element φᵢ
   (`xy_conv_stamp`) and accumulate normal equations (`build_matrix`,
   `build_scprod`)
5. Solve the normal equations via LAPACK Cholesky (`dpotrf` + `dpotrs`) to obtain the cᵢ
   coefficients
6. Sigma-clip bad stamps and refit if necessary
7. Apply the fitted spatially-varying kernel to the full region
   (`spatial_convolve` or `spatial_convolve_fft` if `-fft` flag enabled)
8. Subtract and write output FITS

**Performance (baseline: single-threaded, direct convolution):**

- `spatial_convolve()`: ~60% CPU
- `getPsfCenters()`: ~35% CPU
- `xy_conv_stamp()`: ~18% (subset of above)
- LAPACK solver: <5%

**With optimisations (OpenMP + FFTW3):**

- Region loop parallelisation: 4–8× speedup across cores
- FFT convolution: 3–8× faster on large images (O(N log N) vs. O(k²N))
- Typical total runtime reduction: 8–15× on modern multi-core systems

---

## Building

### CMake (Recommended)

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/hotpants [options]
```

**CMake configuration:**

- Auto-detects CFITSIO, OpenBLAS, FFTW3, and OpenMP
- Optional flags:
    - `-DUSE_FFTW=ON/OFF` — enable FFTW3 convolution (default: ON if found)
    - `-DUSE_OPENMP=ON/OFF` — enable multi-threading (default: ON if found)
    - `-DCMAKE_BUILD_TYPE=Release/Debug` — optimization level

**External dependencies:**

- **CFITSIO** (`-lcfitsio`) — FITS I/O library
- **BLAS/LAPACK** (`-lopenblas` or `-llapack -lblas`) — linear algebra solver
- **FFTW3** (`-lfftw3f`, optional) — FFT-accelerated convolution
- **OpenMP** (`-fopenmp`, optional) — multi-threading support

**Compiler flags:** `-O3 -march=native -funroll-loops -std=c17 -Wall`

Use `-march=native` for SIMD optimisations on the build host; omit for portable
binaries.

### Python (In Development)

```sh
uv sync --dev       # Install dependencies
uv build            # Build wheel and sdist
uv run pytest       # Run tests
```

---

## Future Development

### Python API & Documentation

**Completed:**
- ✓ ctypes bindings with full numpy support
- ✓ Configuration dataclasses (KernelConfig, RegionLayout, NoiseThresholds, KernelSolution)
- ✓ Parameter validation and type hints
- ✓ All unit tests passing

**Optional future enhancements:**
1. Write user guides (Installation, Command-line Reference, Python API Tutorial)
2. Expand `CONTRIBUTING.md` with build, style, testing, and profiling guidance
3. Full end-to-end integration tests with real astronomical images

### Performance Optimisations — Candidate Areas

**PSF-centre finding acceleration (`getPsfCenters()` → ~35% CPU)**

Options:

1. **SIMD vectorisation** — inner stamp-search loops are data-parallel; use
   `-march=native` for autovectorisation or explicit AVX2/AVX-512 intrinsics
2. **Approximate nearest-neighbour** — replace brute-force all-pairs comparison with
   spatial hashing or k-d trees if stamp count is large
3. **GPU offload** — CUDA/OpenCL could handle stamp comparison kernel, but data transfer
   overhead may dominate for modest image sizes

**Convolution kernel basis optimisation**

Current Gaussian basis fixed (σ = 0.7, 1.5, 3.0 px; polynomial degrees 6, 4, 2). For
specific PSF shapes:

- Data-driven basis (PCA of residual PSFs) might reduce `nCompKer` and FFT count
- Explore adaptive basis selection

**Memory and cache**

- `spatial_convolve_fft()` allocates large temporary buffers; consider prefetching or
  tiling strategies
- Monitor cache misses with `perf stat` or `likwid` on representative workloads

**Polynomial coefficient evaluation**

`make_kernel()` evaluates low-order polynomial at every pixel. Current inline; consider:

1. Horner's method (already used) vs. precomputed lookup tables + trilinear
   interpolation
2. Vectorisation of polynomial loop if compiler cannot autovectorise

**Mask propagation in FFT path**

Mask propagation (`makeNoiseImage4()`) currently uses direct convolution even when pixel
values use FFT. Could harmonize to FFT path, but requires careful handling of
integer/binary mask semantics.


---

## Code Quality Standards

### Readability

Core improvements in place:
- Global variable naming legend (`globals.h`): `t*`, `i*`, `m*`, `o*` prefixes
- Named constants in `defaults.h` (histogram, RNG, iteration caps)
- Doxygen `@brief` + `@details` blocks on core functions
- Descriptive loop variables (`pixelX`, `stampCenterX`, `substampIndexX`, etc.)
- Algorithm documentation with Alard & Lupton (1998) references

Future enhancements (optional):
- Single-letter loop variables in nested contexts → descriptive names
- Extract large functions (>150 lines) into testable sub-units
- Naming convention standardization (legacy camelCase → snake_case)

### Contributor Checklist

When writing or modifying code in HOTPANTS, ensure:

- ✓ Global variables follow the `t*/i*/m*/o*` naming convention documented in globals.h
- ✓ Loop variables in nesting depth > 2 use descriptive names (not i, j, k)
- ✓ All magic numbers are #define constants in defaults.h with explanatory comments
- ✓ All functions longer than 50 lines have Doxygen `@brief` and `@details` blocks
- ✓ All algorithms reference Alard & Lupton (1998) or relevant paper (with equations)
- ✓ Non-obvious operations (parity checks, bit-twiddling) have explanatory comments
- ✓ Memory allocation uses `xmalloc()`, `xcalloc()`, or contiguous allocators
- ✓ Output uses `LOG_PROGRESS()`, `LOG_DEBUG()`, `LOG_ERROR()` macros (not direct fprintf)
- ✓ Functions decomposed if nesting depth >3 or cyclomatic complexity >5
- ✓ Code follows conventions: snake_case for new functions, camelCase for legacy/types
- ✓ Performance-critical sections commented with expected CPU share / bottleneck rationale

### Why Readability Matters

1. **Python API Documentation:** C core must be self-documenting for Doxygen + Sphinx;
   vague names and hidden algorithms hinder API docs
2. **Faster Onboarding:** Explicit naming and documented algorithms reduce contributor
   ramp-up time
3. **Debugging & Optimization:** Clear algorithm docs make it easier to profile,
   identify bottlenecks, and verify correctness
4. **Long-term Maintainability:** Well-named code and comprehensive comments reduce
   technical debt and prevent subtle bugs

---

## Development Workflows

### Testing

**Unit tests:**

```sh
uv run pytest tests/test_api.py -v
```

**Integration tests:**

```sh
uv run pytest tests/test_api_integration.py -v
```

**Regression tests:**

```sh
uv run pytest tests/test_regression.py -v
```

**All tests:**

```sh
uv run pytest tests/ -v
```

### Building & Installing

**CMake (C core only):**

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/hotpants [options]
```

**Python package (in development):**

```sh
uv sync --dev         # Install dev dependencies
uv build              # Build wheel and sdist
uv pip install dist/hotpants-*.whl  # Install locally
```

### Documentation

**Generate Sphinx docs (local):**

```sh
cd docs
make html
open _build/html/index.html
```

### Performance Profiling

**Single-threaded baseline:**

```sh
./build/hotpants -iimage.fits -ttemplate.fits -ooutput.fits -n 1
```

**Profile with gprof:**

```sh
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build
OMP_NUM_THREADS=1 ./build/hotpants [options]
gprof ./build/hotpants gmon.out | head -50
```

**Profile with perf:**

```sh
perf record -g ./build/hotpants [options]
perf report
```

### Code Style & Conventions

**Naming:**

- Global variables: `t*`, `i*`, `m*`, `o*` prefixes (documented in globals.h)
- Functions: snake_case for new code, camelCase acceptable for legacy
- Loop variables: descriptive names (pixelX, stampCenterX, etc.) for nesting depth > 2
- Constants: UPPERCASE with #define in defaults.h

**Comments:**

- Doxygen blocks (`@brief`, `@details`, `@param`, `@return`) for functions >50 lines
- Inline comments only for non-obvious logic (algorithms, bit-twiddling, performance
  tricks)
- Reference Alard & Lupton (1998) with equation numbers where relevant

**Compiler flags:**

- `-O3 -march=native -funroll-loops -std=c17 -Wall`
- Use `-march=native` for native binaries; omit for portability
- Enable ASan/UBSan in debug builds for memory safety

---

## Improvements Wishlist

_This section has been added by the user_ - as these are completed, strike them through
or expand them to their own sections as in-progress tasks are completed.

* Implementation of the [Bramich (2008)](https://arxiv.org/abs/0802.1273) delta function
  basis as a selectable option. This directly solves for the kernel, and so generally
  provides more flexibility. This will likely need regularisation to avoid kernel noise,
  via the Laplacian to penalise excessive curvature, or other methods to enforce
  smoothness. This too should ideally be spatially-varying across the frame.
* Implementation of the "afterburner decorrelation"
  from [LSST-DM-TN021](https://dmtn-021.lsst.io/) - this post-convolves the difference
  image with a whitening kernel to remove correlated noise. Some thought is needed on
  how to handle this in a spatially-varying way, since the whitening kernel is not a
  linear combination of basis components.
* Big one: using thin plate splines (or B-splines) to model spatial variation and
  background instead of polynomials. This coupled with setting `-nrx` and`-nry` to 1
  means that entire images can be processed in one go, rather than tiled which can
  introduce step discontinuities in the image. This really shifts the needle for
  widefield instruments
* `SubtractionResult` dataclass which provides Python access to all `hotpants` output
  layers.

---

## Notes for AI Assistants

When working on HOTPANTS:

1. **Understand the algorithm first.** Read Alard & Lupton (1998), Section 2 before
   modifying kernel fitting code.
2. **Prioritize readability.** Clear variable names and algorithm documentation are more
   valuable than micro-optimizations.
3. **Test before optimizing.** Establish baseline performance with gprof/perf before
   attempting speedups.
4. **Respect the modernisation plan.** Python API, Sphinx docs, and Phase 2 readability
   improvements are in progress; coordinate efforts.
5. **Maintain backward compatibility.** CLI and FITS I/O must continue to work; new
   Python API should wrap, not replace, existing C code.
6. **Check the status section.** This file documents what's been completed; avoid
   duplicating work.
7. **Document with Doxygen.** New functions should have `@brief` + `@details` blocks;
   reference papers and equations.
8. **Use named constants.** No magic numbers in code; define #define in defaults.h if
   you need a constant.
9. **Profile on real data.** Optimizations should be measured against representative
   astronomical images.
10. **Communicate with contributors.** If extending the Python API or changing the build
    system, open a discussion in CONTRIBUTING.md or a PR comment.

---

## References

- **Alard & Lupton (1998):** "A Method for Optimal Image Subtraction," ApJ 503:
    325. https://iopscience.iop.org/article/10.1086/305984

    - Section 2: Kernel parameterisation and fitting algorithm
    - Section 3: Convolution and difference image computation
- **CFITSIO:** https://heasarc.gsfc.nasa.gov/fitsio/
- **FFTW3:** https://www.fftw.org/
- **OpenBLAS:** https://www.openblas.net/
- **LAPACK:** https://www.netlib.org/lapack/
- **OpenMP:** https://www.openmp.org/
