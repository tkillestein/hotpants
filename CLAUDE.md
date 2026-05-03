# HOTPANTS — Development Guide

**High Order Transform of PSF ANd Template Subtraction**

HOTPANTS produces difference images between an astronomical science image and a
template by fitting and applying a spatially-varying convolution kernel that
matches the two PSFs. The core algorithm is from Alard & Lupton (1998), ApJ
503:325 (https://iopscience.iop.org/article/10.1086/305984).

License: MIT (Andy Becker, 2013). See `LICENSE`.

---

## Project Status (May 2026)

**Core algorithm & performance:**
- ✓ CMake build system with CFITSIO, OpenBLAS, FFTW3, OpenMP detection
- ✓ LAPACK Cholesky solver (replaced Numerical Recipes LU decomposition)
- ✓ OpenMP parallelism (region loop, FFT inner loops, dynamic scheduling)
- ✓ FFT-accelerated convolution via FFTW3 (O(N log N), 3–8× speedup)

**Code quality & readability:**
- ✓ Doxygen documentation comments added to core functions (`alard.c`, `functions.c`)
- ✓ Global variable naming legend documented in `globals.h` (t*, i*, m*, o* prefixes)
- ✓ Magic numbers replaced with named #define constants in `defaults.h`
- ✓ Histogram algorithm parameters documented (HISTOGRAM_SAMPLE_SIZE, etc.)
- ✓ Loop variables refactored with descriptive names in kernel-fitting code

**Python API & documentation:**
- ✓ Python API header (`hotpants_api.h`) with function signatures and data structures
- ✓ Test structure established (`test_api.py`, `test_api_integration.py`)
- ✓ `pyproject.toml` configured for Python ≥ 3.11, numpy/pytest dependencies
- ✓ Sphinx + Breathe setup for API documentation

**In progress:**
- ⧉ Python CFFI bindings (wrapper code)
- ⧉ End-to-end integration tests
- ⧉ CONTRIBUTING.md guide

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

| Function | File | Purpose | Status |
|---|---|---|---|
| `main()` | `main.c` | Orchestrates the full pipeline | ✓ |
| `fitKernel()` | `alard.c` | Least-squares kernel fit per region | ✓ Documented |
| `build_matrix()` / `build_scprod()` | `alard.c` | Accumulate normal equations | ✓ Documented |
| `dpotrf`/`dpotrs` (LAPACK) | `alard.c` | Cholesky solve of normal equations | ✓ Replaced LU |
| `spatial_convolve()` | `alard.c` | Direct convolution (primary bottleneck, ~60% CPU) | ✓ Documented |
| `spatial_convolve_fft()` | `alard.c` | FFT-accelerated convolution (3–8× faster) | ✓ |
| `make_kernel()` | `alard.c` | Evaluate kernel at image position from polynomials | ✓ Documented |
| `xy_conv_stamp()` | `alard.c` | Separable 2D convolution with Gaussian basis | ✓ Documented |
| `getKernelVec()` | `alard.c` | Initialise multi-Gaussian kernel basis vectors | ✓ Documented |
| `buildStamps()` | `functions.c` | Create stamp grid and compute statistics | ✓ In API |
| `getPsfCenters()` | `functions.c` | Locate bright star centres (secondary bottleneck, ~35% CPU) | ✓ In API |
| `getStampStats3()` | `functions.c` | Sigma-clipped statistics and histogram FWHM | ✓ In API, Documented |
| `makeNoiseImage4()` | `functions.c` | Propagate noise through convolution | ✓ Documented |
| `vargs()` | `vargs.c` | Parse all command-line options | ✓ |

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
3. In each stamp, locate bright PSF centres (`getPsfCenters`) — identifies ~3 substamps per stamp
4. For each substamp, convolve I (or T) with every basis element φᵢ
   (`xy_conv_stamp`) and accumulate normal equations (`build_matrix`,
   `build_scprod`)
5. Solve the normal equations via LAPACK Cholesky (`dpotrf` + `dpotrs`) to obtain the cᵢ coefficients
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

Use `-march=native` for SIMD optimisations on the build host; omit for portable binaries.

### Legacy Makefile

A legacy Makefile is maintained for backward compatibility, but CMake is preferred.

### Python (In Development)

```sh
uv sync --dev       # Install dependencies
uv build            # Build wheel and sdist
uv run pytest       # Run tests
```

---

## Completed Work

### ✓ Linear Algebra — LAPACK Cholesky & CMake

- Replaced `ludcmp()` and `lubksb()` (Numerical Recipes, pre-2000) with LAPACK Cholesky (`dpotrf` + `dpotrs`)
- Normal equations matrix from `build_matrix()` is symmetric positive-definite → Cholesky is optimal (2× fewer operations than general LU)
- Linear system solve now uses optimised library code; typically <5% of runtime
- CMake build system cleanly detects and links CFITSIO, OpenBLAS, FFTW3, and OpenMP
- Both CMake and legacy Makefile available for compatibility

### ✓ Parallelism — OpenMP

- Region loop in `main.c` parallelised: `#pragma omp parallel for schedule(dynamic)`
  - Each region processed independently; dynamic scheduling balances load
- `spatial_convolve_fft()` inner pixel loop parallelised with `collapse(2)` (2D grid)
- FFTW plan creation serialised with `#pragma omp critical(fftw_plan)` (not thread-safe)
- CFITSIO I/O protected with `#pragma omp critical(cfitsio)` (not thread-safe)
- Respects `OMP_NUM_THREADS` at runtime; typical 4–8× thread scaling observed on representative workloads

### ✓ FFT Convolution — FFTW3 Acceleration

- New `spatial_convolve_fft()` function implements O(N log N) convolution via FFTW3
- Restructured kernel as K = Σᵢ cᵢ(x,y)·φᵢ; each basis convolution I⊗φᵢ computed once via FFT, then summed
- FFT sizes optimised to highly-composite numbers (powers of 2, 3, 5) via `next_fftw_size()` for speed
- Plans created once per region, reused across all basis elements
- Mask propagation and variance still computed via direct convolution (unchanged, for correctness)
- Switchable at runtime: `-fft` flag enables FFTW path; default is direct convolution for backward compatibility

### ✓ Readability Improvements — Phase 1

#### Global Variable Naming Legend (`globals.h:34–60`)

Added comprehensive legend documenting all single-letter prefixes:
- `t*` = template image data (e.g., `tUThresh`, `tGain`, `tRdnoise`)
- `i*` = science (input) image (e.g., `iUThresh`, `iGain`)
- `m*` = mask and masking data (e.g., `mRData`)
- `o*` = output image (e.g., `outim`)
- `e*` = ephemeris/environment (reserved)

#### Magic Numbers → Named Constants (`defaults.h:77–100`)

Replaced hardcoded values with documented #define constants:
- **Histogram parameters:** `HISTOGRAM_SAMPLE_SIZE (100)`, `HISTOGRAM_NUM_BINS (256)`, `HISTOGRAM_LOWER_FRAC (0.5)`, `HISTOGRAM_UPPER_FRAC (0.9)`, `HISTOGRAM_PEAK_WIDTH_PCT (0.1)`, `HISTOGRAM_NOISE_HALF_PCT (0.25)`
- **RNG:** `RNG_SEED_MAGIC (-666)` — Numerical Recipes convention; triggers re-seeding of Park & Miller generator
- **Iteration caps:** `MAX_SIGMA_CLIP_RETRIES (5)` — prevents divergence in histogram bin-width refinement

#### Doxygen Documentation Blocks

Added `@brief` and `@details` blocks to core functions in `alard.c` and `functions.c`:
- `fillStamp()` — explains polynomial basis construction and substamp processing
- `fitKernel()` — documents normal equation accumulation and LAPACK solve
- `spatial_convolve()`/`spatial_convolve_fft()` — kernel evaluation and convolution strategy
- `getStampStats3()` — histogram algorithm and sigma-clipping procedure
- References to **Alard & Lupton (1998)** with equation numbers where relevant

#### Loop Variable Refactoring

Renamed cryptic single-letter loop variables with descriptive names:
- `pixelX`, `pixelY` — image pixel coordinates
- `stampCenterX`, `stampCenterY` — stamp position in region
- `substampIndexX`, `substampIndexY` — substamp grid index
- `vectorCompIdx`, `gaussianCompIdx` — vector/Gaussian component indices
- `polyBasisX`, `polyBasisY` — polynomial basis evaluation
- `bgDegX`, `bgDegY` — background polynomial degrees

#### Algorithm Documentation

Added comprehensive comments explaining:
- **Separable convolution:** "Gaussian-polynomial filters are separable: G(x,y)·P(x,y) = [G(x)·P_x(x)] × [G(y)·P_y(y)]. Reduces 2D O(k²n²) to two 1D passes O(kn²), critical for performance."
- **Renormalization flag:** Explains when/why half-Gaussian overlap is handled
- **Parity checks:** Clarifies odd/even detection logic in polynomial basis construction
- **Polynomial basis formula:** Documents triangular basis structure: (order+1)×(order+2)/2 terms

### ✓ Python API — Header & Test Structure

#### Public C API (`hotpants_api.h`)

Designed minimal, Python-friendly interface:
- **Data structures:** `stamp_struct` with clear field documentation
- **Core functions:** `buildStamps()`, `getPsfCenters()`, `fitKernel()`, `spatial_convolve()`, `getStampStats3()`, etc.
- **Global configuration:** `hwKernel`, `kerOrder`, `bgOrder`, `nCompKer`, `nComp`, `verbose`, `nThread`, etc.
- Each function includes Doxygen `@brief`, parameter descriptions, and reference to Alard & Lupton (1998)
- All fields and functions documented for automatic API doc generation

#### Test Structure (`tests/`)

- `conftest.py` — pytest fixtures and configuration
- `test_api.py` — unit tests for `KernelConfig`, `RegionLayout`, `NoiseThresholds`, `KernelSolution` dataclasses
  - Configuration validation, default values, parameter ranges
- `test_api_integration.py` — end-to-end integration tests (TBD)
- `test_regression.py` — regression tests against known outputs (TBD)

#### Python Configuration (`pyproject.toml`)

- Requires Python ≥ 3.11
- Dependencies: `numpy >= 1.26`
- Dev dependencies: `pytest >= 8.0`, `numpy`, `astropy`, `scipy`
- Docs dependencies: `sphinx`, `sphinx-rtd-theme`, `breathe`
- Build backend: `hatchling` (compatible with `uv build`)

### ✓ Documentation Setup — Sphinx & Breathe

- **Sphinx configuration** (`docs/conf.py`) with `breathe` extension for Doxygen integration
- **RTD theme** for modern, responsive documentation
- **Source structure:** `docs/api/`, `docs/guides/`, `docs/index.rst`
- **Doxygen output** auto-included in Sphinx build

---

## Modernisation Goals — Ongoing

This section records the intended development direction. Follow these
preferences in all future work on the repository.

### Python Bindings — CFFI (In Progress)

**Current status:** API header (`hotpants_api.h`) complete; Python test structure in place.

**Next steps:**
1. Write CFFI build wrapper (see `pyproject.toml` for hatchling integration)
2. Implement high-level Python module (`hotpants/__init__.py`) with:
   - `KernelConfig` dataclass (kernel_half_width, kernel_order, bg_order, etc.)
   - `fit_kernel()` function accepting numpy arrays
   - `spatial_convolve()` function
   - Error handling and validation
3. Add numpy array marshalling (C `float*` ↔ `numpy.ndarray` conversion)
4. Write integration tests covering full workflow (FITS load → fit → convolve → validate)

**Design principles:**
- Accept and return **numpy arrays** (no FITS I/O in C extension)
- Callers use `astropy.io.fits` for FITS I/O
- Keep Python layer thin; heavy lifting in optimised C
- Parameter validation in Python; error handling from C via return codes

### Documentation — Sphinx/Breathe (In Progress)

**Current status:** Sphinx + Breathe configured; Doxygen comments in place.

**Next steps:**
1. Generate Doxygen XML output (run `doxygen` in build)
2. Write user guides:
   - "Installation & Quick Start" (pip/uv commands)
   - "Command-line Reference" (copy/expand from README.md)
   - "Python API Tutorial" (image arrays → difference image workflow)
   - "Algorithm Overview" (Alard & Lupton theory + implementation)
3. Write `CONTRIBUTING.md` covering:
   - Build instructions (CMake setup, dependencies)
   - Code style (snake_case for new functions, Doxygen comments)
   - Testing expectations (unit tests for new features, pytest before pushing)
   - Performance profiling (gprof, perf, cache analysis)

### Performance Optimisations — Future

**PSF-centre finding acceleration (`getPsfCenters()` → ~35% CPU)**

Options:
1. **SIMD vectorisation** — inner stamp-search loops are data-parallel; use `-march=native` for autovectorisation or explicit AVX2/AVX-512 intrinsics
2. **Approximate nearest-neighbour** — replace brute-force all-pairs comparison with spatial hashing or k-d trees if stamp count is large
3. **GPU offload** — CUDA/OpenCL could handle stamp comparison kernel, but data transfer overhead may dominate for modest image sizes

**Convolution kernel basis optimisation**

Current Gaussian basis fixed (σ = 0.7, 1.5, 3.0 px; polynomial degrees 6, 4, 2). For specific PSF shapes:
- Data-driven basis (PCA of residual PSFs) might reduce `nCompKer` and FFT count
- Explore adaptive basis selection

**Memory and cache**

- `spatial_convolve_fft()` allocates large temporary buffers; consider prefetching or tiling strategies
- Monitor cache misses with `perf stat` or `likwid` on representative workloads

**Polynomial coefficient evaluation**

`make_kernel()` evaluates low-order polynomial at every pixel. Current inline; consider:
1. Horner's method (already used) vs. precomputed lookup tables + trilinear interpolation
2. Vectorisation of polynomial loop if compiler cannot autovectorise

**Mask propagation in FFT path**

Mask propagation (`makeNoiseImage4()`) currently uses direct convolution even when pixel values use FFT. Could harmonize to FFT path, but requires careful handling of integer/binary mask semantics.

### Code Quality — Ongoing

**Priority readability improvements (remaining items from initial audit):**

See "Readability Improvements — Phase 2" section below.

---

## Readability Improvements — Status & Roadmap

This section tracks code clarity improvements. Phase 1 (magic numbers, global naming, Doxygen comments) is **COMPLETE**. Remaining items are lower-priority but important for long-term maintainability.

### Phase 1 — Completed ✓

- ✓ Global variable naming legend (globals.h:34–60)
- ✓ Magic numbers → #define constants (defaults.h:77–100)
- ✓ Doxygen `@brief` + `@details` blocks on core functions
- ✓ Loop variable refactoring (pixelX, stampCenterX, substampIndexX, etc.)
- ✓ Algorithm documentation (separable convolution, renormalization, polynomial basis)

### Phase 2 — Medium Priority (2–6 hours per item)

| Item | Status | Effort | Notes |
|---|---|---|---|
| Single-letter loop variables in nested contexts (i, j, k, l, m, n) | Pending | 2–3 hrs | Apply to `alard.c`, `functions.c` pervasively; use `pixel_x`, `bin_idx`, `vector_idx`, `stamp_idx` based on context |
| Cryptic variable names (`q`, `ncomp1`, `ncomp2`, `ivecbg`) | Pending | 1–2 hrs | Replace with full names: `scalar_product`, `num_kernel_components_low`, `num_kernel_components_high`, `background_vector_index` |
| Separate large functions | Pending | 4–6 hrs | Extract `extract_stamp_matrix()`, `solve_kernel_linear_system()`, `compute_sigma_clipped_mean()`, `estimate_histogram_fwhm()` from monolithic functions |
| Histogram FWHM algorithm extraction | Pending | 2 hrs | Extract `getStampStats3()` logic into helper functions with clear single responsibility |
| Separable convolution documentation | Pending | 30 min | Expand comment on why separability is chosen (O(k²n²) → O(kn²) reduction) |

### Phase 3 — Lower Priority (Nice-to-Have)

| Item | Status | Effort | Notes |
|---|---|---|---|
| Dead code removal | Pending | 1–2 hrs | Remove `/* DD fprintf(...) */` fragments and obsolete variable definitions |
| Consistent logging macro | Pending | 2–3 hrs | Standardise `verbose >= 1/2` checks with centralized `VERBOSE_LOG(level, fmt, ...)` macro |
| Naming convention standardization | Ongoing | 8–10 hrs | Migrate legacy camelCase to snake_case; document in CONTRIBUTING.md |
| Full function decomposition | Pending | 6–8 hrs | Break oversized functions (>150 lines) into testable sub-units |

### Contributor Checklist

When writing or modifying code in HOTPANTS, ensure:

- ✓ Global variables follow the `t*/i*/m*/o*` naming convention documented in globals.h
- ✓ Loop variables in nesting depth > 2 use descriptive names (not i, j, k)
- ✓ All magic numbers are #define constants in defaults.h with explanatory comments
- ✓ All functions longer than 50 lines have Doxygen `@brief` and `@details` blocks
- ✓ All algorithms reference Alard & Lupton (1998) or relevant paper (with equation numbers)
- ✓ Non-obvious operations (parity checks, bit-twiddling) have explanatory comments
- ✓ Functions decomposed if nesting depth >3 or cyclomatic complexity >5
- ✓ Dead code removed; debug output uses consistent verbose-level checks
- ✓ Code follows conventions: snake_case for new functions, camelCase reserved for legacy/types
- ✓ Performance-critical sections commented with expected CPU share / bottleneck rationale

### Why Readability Matters

1. **Python API Documentation:** C core must be self-documenting for Doxygen + Sphinx; vague names and hidden algorithms hinder API docs
2. **Faster Onboarding:** Explicit naming and documented algorithms reduce contributor ramp-up time
3. **Debugging & Optimization:** Clear algorithm docs make it easier to profile, identify bottlenecks, and verify correctness
4. **Long-term Maintainability:** Well-named code and comprehensive comments reduce technical debt and prevent subtle bugs

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
- Inline comments only for non-obvious logic (algorithms, bit-twiddling, performance tricks)
- Reference Alard & Lupton (1998) with equation numbers where relevant

**Compiler flags:**
- `-O3 -march=native -funroll-loops -std=c17 -Wall`
- Use `-march=native` for native binaries; omit for portability
- Enable ASan/UBSan in debug builds for memory safety

---

## Notes for AI Assistants

When working on HOTPANTS:

1. **Understand the algorithm first.** Read Alard & Lupton (1998), Section 2 before modifying kernel fitting code.
2. **Prioritize readability.** Clear variable names and algorithm documentation are more valuable than micro-optimizations.
3. **Test before optimizing.** Establish baseline performance with gprof/perf before attempting speedups.
4. **Respect the modernisation plan.** Python API, Sphinx docs, and Phase 2 readability improvements are in progress; coordinate efforts.
5. **Maintain backward compatibility.** CLI and FITS I/O must continue to work; new Python API should wrap, not replace, existing C code.
6. **Check the status section.** This file documents what's been completed; avoid duplicating work.
7. **Document with Doxygen.** New functions should have `@brief` + `@details` blocks; reference papers and equations.
8. **Use named constants.** No magic numbers in code; define #define in defaults.h if you need a constant.
9. **Profile on real data.** Optimizations should be measured against representative astronomical images.
10. **Communicate with contributors.** If extending the Python API or changing the build system, open a discussion in CONTRIBUTING.md or a PR comment.

---

## References

- **Alard & Lupton (1998):** "A Method for Optimal Image Subtraction," ApJ 503:325. https://iopscience.iop.org/article/10.1086/305984
  - Section 2: Kernel parameterisation and fitting algorithm
  - Section 3: Convolution and difference image computation
- **CFITSIO:** https://heasarc.gsfc.nasa.gov/fitsio/
- **FFTW3:** https://www.fftw.org/
- **OpenBLAS:** https://www.openblas.net/
- **LAPACK:** https://www.netlib.org/lapack/
- **OpenMP:** https://www.openmp.org/
