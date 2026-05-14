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
- ✓ LAPACK Cholesky solver, OpenMP parallelism, FFT-accelerated convolution (mandatory)
- ✓ Tier-2 performance optimizations: cache-aware tiling, SIMD vectorization
- ✓ Doxygen-documented core functions, global naming legend, named constants
- ✓ Logging macros, memory allocation wrappers, contiguous allocation
- ✓ OpenMP parallelization of region and block loops with proper thread-local state

**Python API (ctypes bindings):**
- ✓ C API header (`hotpants_api.h`) with function signatures and data structures
- ✓ High-level Python API: `fit_kernel()`, `spatial_convolve()` with full numpy support
- ✓ Configuration dataclasses: `KernelConfig`, `RegionLayout`, `NoiseThresholds`, `KernelSolution`
- ✓ Pydantic validation with informative error messages
- ✓ All unit tests passing (33/33), segfaults fixed, parameter validation working
- ✓ Logging via loguru with integration to C library output

**Documentation & build:**
- ✓ `pyproject.toml` configured for Python ≥ 3.11 (excludes 3.14 due to segfault)
- ✓ Sphinx + Breathe setup, Doxygen integration
- ✓ Test suite: unit tests, integration tests, regression tests
- ✓ Performance profiling tools: `profile_hotpants.py`, `analyze_bottlenecks.py`
- ✓ `BOTTLENECK_ANALYSIS.md`: comprehensive profiling report with optimization roadmap

---

## Repository Structure

```
hotpants/
├── src/
│   ├── main.c              # Pipeline orchestration: FITS I/O, region loop, output assembly
│   ├── alard.c             # Core algorithm: kernel fitting, spatial convolution, LAPACK solve
│   │                       # (2470 lines: contains multiple code paths for optimization)
│   ├── functions.c         # Stamps, PSF-centre finding, statistics, masking (2343 lines)
│   ├── hotpants_wrapper.c  # C wrapper for Python API; global state setup (528 lines)
│   ├── vargs.c             # CLI argument parsing (~50 options, 754 lines)
│   ├── maskim.c            # Standalone utility: apply mask to image
│   ├── globals.c           # Global variable definitions (11 lines)
│   ├── allocate.c          # Memory allocation wrappers (147 lines)
│   ├── defaults.h          # Compile-time parameter defaults (with named constants)
│   ├── globals.h           # Global variable declarations + prefix legend
│   ├── allocate.h          # Allocation wrapper prototypes
│   ├── functions.h         # Function prototypes and CFITSIO includes
│   └── hotpants_api.h      # Public C API for Python bindings (Doxygen-documented)
├── src/hotpants/           # Python package
│   ├── __init__.py         # Public API exports
│   ├── api.py              # High-level API (fit_kernel, spatial_convolve, config dataclasses)
│   ├── _core.py            # Low-level ctypes wrappers, array conversion, global state
│   └── _hotpants_ffi.py    # FFI bindings to compiled C library
├── tests/
│   ├── conftest.py         # pytest configuration and test fixtures
│   ├── test_api.py         # Unit tests for Python API (config validation, error handling)
│   ├── test_api_integration.py  # End-to-end integration tests with synthetic data
│   ├── test_regression.py  # Regression tests against known outputs
│   └── __init__.py
├── docs/
│   ├── conf.py             # Sphinx configuration with Breathe (Doxygen integration)
│   ├── index.rst           # Documentation root
│   ├── api/                # API reference (generated from Doxygen + docstrings)
│   ├── guides/             # User guides and tutorials
│   ├── _static/            # Static assets (CSS, images)
│   └── _templates/         # Sphinx templates
├── CMakeLists.txt          # CMake build configuration (FFTW3 mandatory)
├── Doxyfile                # Doxygen configuration for C API documentation
├── pyproject.toml          # Python project metadata (uv-managed)
├── README.md               # User-facing option reference and quick start
├── CLAUDE.md               # This file: development guide for AI assistants
├── BOTTLENECK_ANALYSIS.md  # Profiling report and acceleration roadmap
├── CONTRIBUTING.md         # Contributor guidelines
├── NOTES.md                # Version history and changelog
├── profile_hotpants.py     # Profiling utility for performance analysis
├── analyze_bottlenecks.py  # Bottleneck identification and reporting
└── profile_gprof.py        # gprof output parsing and analysis
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

### System Dependencies

Before building, install required system libraries:

**Debian/Ubuntu:**
```sh
sudo apt-get install libcfitsio-dev libfftw3-dev liblapack-dev libblas-dev libopenblas-dev liblapacke-dev
```

**macOS (using Homebrew):**
```sh
brew install cfitsio fftw openblas lapack
```

**RHEL/CentOS/Fedora:**
```sh
sudo dnf install cfitsio-devel fftw-devel lapack-devel openblas-devel
```

**Required Libraries:**
- **CFITSIO** (`libcfitsio-dev`): FITS file I/O
- **FFTW3** (`libfftw3-dev`): Fast Fourier Transform (mandatory for FFT convolution)
- **BLAS** (`libblas-dev` or `libopenblas-dev`): Basic Linear Algebra
- **LAPACK** (`liblapack-dev`): Dense linear algebra solver
- **LAPACKE** (`liblapacke-dev`): C interface for LAPACK

### CMake (Recommended)

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/hotpants [options]
```

**CMake configuration:**

- Auto-detects CFITSIO, OpenBLAS, FFTW3, and OpenMP
- FFTW3 is now required (FFT convolution is always enabled)
- Optional flags:
    - `-DUSE_OPENMP=ON/OFF` — enable multi-threading (default: ON if found)
    - `-DCMAKE_BUILD_TYPE=Release/Debug` — optimization level

**External dependencies (linked automatically):**

- **CFITSIO** — FITS I/O library
- **BLAS/LAPACK** — Linear algebra solvers
- **FFTW3** — FFT-accelerated convolution (now mandatory)
- **OpenMP** — Multi-threading (optional)

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

## Recent Performance Work

### Completed Optimizations (Tier 1 & 2)

As of May 2026, the following improvements have been implemented and tested:

**Tier 1: Mandatory FFT path (∼3–8× speedup)**
- FFT-accelerated convolution is now the primary algorithm (direct convolution removed)
- FFTW3 is mandatory at build time; no `#ifdef USE_FFTW` fallback
- FFT plan creation and caching optimized to minimize overhead
- Comprehensive documentation in `BOTTLENECK_ANALYSIS.md` for troubleshooting

**Tier 2: SIMD vectorization & tiling (∼1.5–2× speedup)**
- Inner kernel-summation loops now use `-march=native` for autovectorization
- Cache-aware tiling in `spatial_convolve_fft()` to reduce L1/L2 misses
- OpenMP parallelization of region and block loops with thread-local state
- All tier-2 changes integrated and validated; marked complete in BOTTLENECK_ANALYSIS.md

**Combined speedup:** ∼5–15× on typical multi-core systems (FFT 3–8×, parallelism 4–8×, tiling 1.5–2×)

### Remaining Performance Opportunities (Tier 3)

See `BOTTLENECK_ANALYSIS.md` for detailed analysis of remaining bottlenecks.

**PSF-centre finding (`getPsfCenters()` → ~35% CPU baseline, reduced by tiling)**

Remaining acceleration paths:

1. **GPU offload** — CUDA/OpenCL for all-pairs PSF center comparison, but data transfer
   overhead likely exceeds benefit for typical image sizes
2. **Approximate nearest-neighbour** — spatial hashing or k-d trees for large stamp counts
   (>1000), but adds algorithmic complexity with marginal gains on typical workloads
3. **Specialized hardware** — FPGA matrix multiplication, unlikely cost-effective

**Algorithm research opportunities (exploratory)**

- **Data-driven Gaussian basis** — PCA on residual PSFs to auto-select optimal basis
- **Bramich (2008) delta-function basis** — direct kernel solving with Laplacian smoothness
- **LSST decorrelation** — post-convolve whitening kernel to remove correlated noise
- **Thin-plate splines** — replace polynomial spatial variation, enable single full-frame region

### Documentation & Testing

**Completed:**
- ✓ Comprehensive Python API with numpy marshalling and error handling
- ✓ Configuration dataclasses with Pydantic validation
- ✓ Unit, integration, and regression test suites
- ✓ BOTTLENECK_ANALYSIS.md with profiling methodology and roadmap
- ✓ Performance profiling scripts: `profile_hotpants.py`, `analyze_bottlenecks.py`

**Optional enhancements:**
- User guides (Installation, Python API Tutorial, Command-line Reference)
- Expand CONTRIBUTING.md with profiling and optimization workflow
- Integration tests with real astronomical images


---

## Technical Debt & Refactoring Opportunities

### Summary of Completed Work (May 2026)

**Major cleanup completed:**
- ✓ PCA code stubs removed (~130 lines)
- ✓ `build_matrix_new()` removed (114 lines) — experimental/unused implementation
- ✓ `calculateAvgNoise()` removed (34 lines) — completely unused function
- ✓ `getPsfCentersORIG()` declaration removed (1 line) — orphaned prototype
- ✓ Custom quicksort replaced with POSIX `qsort()` (~23 line net savings)
- **Total cleanup: ~230 lines of dead/unused code removed or refactored**

All changes verified: 58/58 tests pass, build succeeds, no functionality lost.

---

### C Code: Redundant and Parallel Implementations

**Priority 1: Multiple matrix-building functions (reduce from 3 to 1)**

| Function | File | Lines | Status | Issue |
|----------|------|-------|--------|-------|
| `build_matrix()` | alard.c | ~1165–1730 | ✓ Active | Primary implementation; correct |
| `build_matrix0()` | alard.c | ~638–703 | Legacy | Superseded; kept for reference per comment |
| `build_matrix_new()` | alard.c | ~1027–1164 | ✓ **REMOVED** | Experimental version; never called — deleted May 2026 |

**Status:** `build_matrix_new()` has been removed (114 lines). Consider consolidating `build_matrix0()` 
and `build_scprod0()` if the per-stamp code path is no longer needed.

---

**Priority 1: PCA alternative code paths** ✓ **COMPLETED**

| Component | Implementations | Status |
|-----------|-----------------|--------|
| Kernel basis | `kernel_vector_PCA()` | ✓ **REMOVED** — was ~40 lines |
| Stamp convolution | `xy_conv_stamp_PCA()` | ✓ **REMOVED** — was ~50 lines |
| Global variables | `usePCA`, `fwKernelPCA`, `PCA` | ✓ **REMOVED** |
| Fill stamp logic | Dispatcher in `fillStamp()` | ✓ **REMOVED** |

**Status (May 2026):** PCA code stubs have been completely removed. The PCA code was never 
exposed to CLI or Python API, had no tests, and was never initialized from user-facing interfaces. 
Total removed: ~130 lines of dead code.

---

**Priority 2: Float-specific utility duplication** — Partially simplified

| Function | File | Purpose | Status |
|----------|------|---------|--------|
| `fset()` | functions.c | Fill float array | Active; ~6 lines |
| `dfset()` | functions.c | Fill double array | Simplified; ~5 lines (removed pointer) |

**Status:** Functions remain separate for type safety. `dfset()` simplified by removing 
intermediate pointer variable. Full consolidation via macro deferred (minimal savings; type-safe 
separation preferred).

---

**Priority 2: Quick sort implementation** ✓ **COMPLETED**

| Function | File | Status | Note |
|----------|------|--------|------|
| `quick_sort()` | functions.c | ✓ **REPLACED** with qsort() | Was 25 lines; now ~9 lines |
| `quick_sort_1()` | functions.c | ✓ **REMOVED** | Was 59 lines of recursion |
| Helper comparator | functions.c | ✓ **ADDED** | Static indirect comparator: ~4 lines |

**Status (May 2026):** Replaced custom recursive quicksort with POSIX `qsort()` using 
indirect comparator. Net savings: ~23 lines. No performance impact; improves portability 
and reduces maintenance burden.

---

### Python Code: Potential Simplifications

**Priority 2: Parameter duplication in `api.py` (config → global state)**

The high-level API (`fit_kernel()`, `spatial_convolve()`) accepts configuration objects
and translates them to C globals via `global_state()` context manager in `_core.py`.
This works but creates:

1. Dual representation (Python config + C globals)
2. Hidden state management (globals are modified, then restored)
3. Potential race conditions if Python code is multithreaded and calls HOTPANTS multiple times

**Current mitigations:**
- `global_state()` is a context manager ensuring cleanup
- Tests cover single-threaded usage; no multi-threaded tests exist

**Action (exploratory):**
- Consider passing config struct through to C API instead of modifying globals
- Would require changes to C API but would eliminate global state
- Deferred unless thread-safety issues emerge

---

### Repository Structure: Organizational Opportunities

**Priority 3: Separate profiling and analysis tools**

Files `profile_hotpants.py`, `analyze_bottlenecks.py`, `profile_gprof.py` are useful for development
but clutter the root directory. Consider moving to `dev/profiling/` subdirectory.

**Priority 3: BOTTLENECK_ANALYSIS.md maintenance**

This file was generated May 8, 2026. As code changes, it will become stale. Consider:
- Adding a note at the top with last-update date and regeneration instructions
- Automating regeneration as part of CI/CD (gprof on release builds)
- Linking from CONTRIBUTING.md for visibility

---

### Known Limitations & Workarounds

**Python 3.14 segfault**

- Excludes Python 3.14+ from `pyproject.toml` due to segfault in ctypes interaction
- Root cause unknown; likely interaction with numpy/ctypes on new Python release
- Workaround: use Python 3.11–3.13

**Global state thread safety**

- The C library uses global variables for configuration (e.g., `hwKernel`, `kerOrder`)
- OpenMP parallelism is internal to individual function calls, not across calls
- Python API is safe for multithreaded use if each thread serializes HOTPANTS calls (use a mutex)
- If future work adds thread pool support, refactor to pass config struct through C API

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
   modifying kernel fitting code. Check BOTTLENECK_ANALYSIS.md for profiling context.

2. **Prioritize readability over micro-optimizations.** The codebase has already been
   profiled and tier-1/2 optimizations completed. Further speedups must be weighed
   against maintainability.

3. **Check Technical Debt section.** Several code paths (PCA variants, legacy matrix builders)
   are candidates for refactoring. Avoid adding new features that compound duplication.

4. **Test performance impact.** Use `profile_hotpants.py` and `analyze_bottlenecks.py` to
   establish baseline before and after changes. Results should be documented.

5. **Respect backward compatibility.** CLI and FITS I/O must continue to work. New Python
   API features should extend, not replace, existing C signatures.

6. **FFT is mandatory.** Do not introduce code paths that fall back to direct convolution
   or bypass FFTW3. The `.gitignore` excludes FFTW3 library code; ensure system dependencies
   are installed (see "Building" section).

7. **Global state is constrained.** C library uses global variables for config; Python API
   manages them via `global_state()` context manager. If adding new globals, update both
   `_core.py` and `hotpants_wrapper.c`.

8. **Document with Doxygen.** New functions >50 lines need `@brief` + `@details` blocks
   with parameter and return documentation. Reference papers with equation numbers.

9. **Use named constants.** No magic numbers in code; define `#define` in `defaults.h` with
   explanatory comments. Update the global naming legend in `globals.h` if adding new globals.

10. **Profiling workflow:**
    - Baseline: `cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo && cmake --build build`
    - Profile: `OMP_NUM_THREADS=1 python profile_hotpants.py <test_image>`
    - Analyze: `python analyze_bottlenecks.py <gprof_output>`
    - Compare: Document baseline vs. optimized runtime and CPU %; update BOTTLENECK_ANALYSIS.md

11. **Before refactoring legacy code** (e.g., PCA path, matrix builders):
    - Check git log for context on why alternatives were kept
    - Run full test suite (`pytest tests/ -v`) to ensure no regressions
    - Update CLAUDE.md Technical Debt section if status changes

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
