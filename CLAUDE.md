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
- ✓ **Thin Plate Spline (TPS) spatial variation** for kernel and background coefficients
- ✓ **FFTW3 and BLAS multi-threading support** for single-region mode performance (May 2026)
- ✓ **LSST DMTN-021 afterburner decorrelation** with spatially-varying whitening kernel (May 2026)

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
│   ├── decorrelation.c     # DMTN-021 afterburner decorrelation (~750 lines)
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
│   ├── decorrelation.h     # DMTN-021 decorrelation API (Doxygen-documented)
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
spatial interpolation functions:

```
K(x,y) = Σᵢ cᵢ(x,y) · φᵢ
```

where φᵢ are fixed Gaussian profiles (default: 3 Gaussians with σ = 0.7, 1.5,
3.0 pixels and polynomial degrees 6, 4, 2) and cᵢ(x,y) are spatial variation
functions. This parameterization follows **Alard & Lupton (1998), Section 2.2**.

**Spatial variation modes:**

1. **Polynomial mode (default):** cᵢ(x,y) are low-order 2-D polynomials in image position
2. **TPS mode (via `-useTPS 1`):** cᵢ(x,y) and background are Thin Plate Spline (RBF)
   surfaces fitted to stamp-center values, enabling smooth single-region processing

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
sudo apt-get install libcfitsio-dev libfftw3-dev libfftw3-omp liblapack-dev libopenblas-dev liblapacke-dev
```

**macOS (using Homebrew):**
```sh
brew install cfitsio fftw openblas lapack
# Note: fftw on macOS includes OMP support by default
```

**RHEL/CentOS/Fedora:**
```sh
sudo dnf install cfitsio-devel fftw-devel fftw-libs-omp lapack-devel openblas-devel
```

**Required Libraries:**
- **CFITSIO** (`libcfitsio-dev`): FITS file I/O
- **FFTW3** (`libfftw3-dev`): Fast Fourier Transform (mandatory for FFT convolution)
- **FFTW3-OMP** (`libfftw3-omp`): Multi-threaded FFT execution (optional but recommended)
- **OpenBLAS** (`libopenblas-dev`): Multi-threaded linear algebra (preferred over Netlib BLAS)
- **LAPACK** (`liblapack-dev`): Dense linear algebra solver
- **LAPACKE** (`liblapacke-dev`): C interface for LAPACK

**Performance note:** For single-region mode (`-nrx 1 -nry 1`) with TPS, both libfftw3-omp
and libopenblas are essential for 3–8× speedup on multi-core systems.

### CMake (Recommended)

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/hotpants [options]
```

**CMake configuration:**

- Auto-detects CFITSIO, OpenBLAS, FFTW3, FFTW3-OMP, and OpenMP
- FFTW3 is required (FFT convolution is mandatory)
- Prefers OpenBLAS over Netlib BLAS for runtime thread control
- Looks for libfftw3_omp; warns if not available
- Optional flags:
    - `-DUSE_OPENMP=ON/OFF` — enable region-level parallelism (default: ON if found)
    - `-DCMAKE_BUILD_TYPE=Release/Debug` — optimization level

**External dependencies (linked automatically):**

- **CFITSIO** — FITS I/O library
- **FFTW3** — FFT-accelerated convolution (mandatory)
- **FFTW3-OMP** — Multi-threaded FFT execution (optional; improves single-region performance)
- **OpenBLAS** — Multi-threaded BLAS/LAPACK (optional; improves Cholesky solve)
- **OpenMP** — Region-level parallelism (optional; improves multi-region performance)

**Threading recommendations:**

- Single-region mode (`-nrx 1 -nry 1`): Enable both libfftw3_omp and libopenblas for 3–8× speedup
- Multi-region mode: Region-level OpenMP dominates; BLAS/FFTW3 threading provides fallback
- Verify threading is available: `ldd ./build/hotpants | grep -E 'fftw3_omp|openblas'`

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

**Tier 2b: FFTW3 and BLAS/LAPACK multi-threading (∼1.5–4× speedup, single-region mode)**

Completed (May 2026):
- ✓ **FFTW3 multi-threaded FFT execution** — `fftw_init_threads()` + `fftw_plan_with_nthreads()`
  - Requires linking `libfftw3_omp`; automatically detects OpenMP thread count
  - Single-region mode uses `FFTW_MEASURE` for 5–10% better FFT execution
  - Multi-region mode uses `FFTW_ESTIMATE` for fast planning
  - Plans cached by (nx, ny) to avoid replanning; cache persists for program lifetime
  - Estimated speedup: 1.5–2× on multi-core FFT-heavy workloads

- ✓ **BLAS/LAPACK multi-threading** — Detects OpenBLAS or MKL at runtime
  - `init_threading()` called after OpenMP thread pool setup
  - Sets BLAS/LAPACK to `omp_get_max_threads()` for optimal parallelism
  - Cholesky solve (`dpotrf` + `dpotrs`) parallelizes well for large matrices (n > 100)
  - Estimated speedup: 2–4× for typical kernel solutions (nCompKer × nS ≈ 600×600 matrix)

**Single-region mode total speedup** (combining region-level OpenMP + FFTW/BLAS threading):
- With 4 cores: ~3–8× speedup over single-threaded baseline
- Breakdown: FFT 1.5–2×, BLAS 2–4×, combined avoids oversubscription via thread pool control

**Multi-region mode unchanged:** Region-level OpenMP parallelism dominates (4–8×); FFTW/BLAS
threads per region reduced to 1 to avoid oversubscription. Plan caching still provides ~1–2%
benefit if regions have uniform sizes.

### Remaining Performance Opportunities (Tier 3)

See `BOTTLENECK_ANALYSIS.md` for detailed analysis of remaining bottlenecks.

**PSF-centre finding (`getPsfCenters()` → ~35% CPU baseline, reduced by tiling)**

Remaining acceleration paths:

1. **GPU offload** — CUDA/OpenCL for all-pairs PSF center comparison, but data transfer
   overhead likely exceeds benefit for typical image sizes
2. **Approximate nearest-neighbour** — spatial hashing or k-d trees for large stamp counts
   (>1000), but adds algorithmic complexity with marginal gains on typical workloads
3. **Specialized hardware** — FPGA matrix multiplication, unlikely cost-effective

**Algorithm improvements (completed & exploratory)**

**Completed (May 2026):**
- ✓ **Thin Plate Spline (TPS) spatial variation** — Radial basis function (RBF) surfaces
  replace polynomial spatial variation for both kernel coefficients and background.
  Enables smooth full-frame processing with `-nrx 1 -nry 1` (no tiling artifacts).
  Usage: `-useTPS 1 -tpsSmoothing <value>` (default smoothing: 1e-6)

**Remaining research opportunities (exploratory):**

- **Data-driven Gaussian basis** — PCA on residual PSFs to auto-select optimal basis
- **Bramich (2008) delta-function basis** — direct kernel solving with Laplacian smoothness
- **LSST decorrelation** — post-convolve whitening kernel to remove correlated noise

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

**Major algorithmic improvements:**
- ✓ **Thin Plate Spline (TPS) implementation** — complete RBF-based spatial variation system
  for kernel and background (Phases 1–6, ~800 lines of new code)
  - Core RBF infrastructure: `tps_kernel()`, `tps_assemble_matrix()`, `tps_fit_coefficients()`, `tps_evaluate()`
  - Kernel TPS: `tps_fit_kernel()`, `make_kernel_tps()`
  - Background TPS: `tps_fit_background()`, `get_background_tps()`
  - Python API: `KernelConfig.use_tps`, `KernelConfig.tps_smoothing`
  - CLI: `-useTPS`, `-tpsSmoothing` flags
  - All 84 unit/integration/regression tests passing

**Major cleanup completed:**
- ✓ PCA code stubs removed (~130 lines)
- ✓ `build_matrix_new()` removed (114 lines) — experimental/unused implementation
- ✓ `calculateAvgNoise()` removed (34 lines) — completely unused function
- ✓ `getPsfCentersORIG()` declaration removed (1 line) — orphaned prototype
- ✓ Custom quicksort replaced with POSIX `qsort()` (~23 line net savings)
- **Total cleanup: ~230 lines of dead/unused code removed or refactored**

All changes verified: 84/84 tests pass, build succeeds, no functionality lost.

---

### TPS Implementation Architecture (May 2026)

**Overview:**

The Thin Plate Spline (TPS) spatial variation system provides smooth RBF-based
interpolation of kernel coefficients and background, replacing polynomial spatial
variation when enabled.

**Core Components:**

1. **RBF Infrastructure** (alard.c: lines 35–240)
   - `tps_kernel()`: φ(r) = r² log(r) RBF kernel
   - `tps_assemble_matrix()`: Build augmented RBF+polynomial matrix
   - `tps_fit_coefficients()`: LAPACK LU solve for RBF weights + polynomial trend
   - `tps_evaluate()`: Evaluate RBF surface at query point

2. **Kernel TPS** (alard.c: lines 280–384)
   - `tps_fit_kernel()`: Fit RBF surface to each kernel component post-LAPACK solve
   - `make_kernel_tps()`: Evaluate spatially-varying kernel using TPS at pixel (xi, yi)
   - `make_kernel_dispatch()`: Router between TPS and polynomial evaluation based on `useTPS`

3. **Background TPS** (alard.c: lines 386–480)
   - `tps_fit_background()`: Fit RBF surface to background polynomial post-LAPACK solve
   - `get_background_tps()`: Evaluate background value using TPS
   - Offset functions: `kernelSol_offset_tps_bg_weights()`, `kernelSol_offset_tps_bg_poly()`

4. **Integration** (alard.c: lines 815–890)
   - `fitKernel()`: Calls `tps_fit_kernel()` and `tps_fit_background()` when `useTPS=1`
   - `get_background()`: Dispatches to `get_background_tps()` if `useTPS=1`

**kernelSol Layout (TPS Mode):**

```
Index Range                           | Content
[0..poly_size-1]                      | Polynomial coefficients (for compat)
[poly_size..poly_size+nCompKer*nS-1] | Kernel RBF weights (nCompKer components)
[..+3*nCompKer-1]                     | Kernel polynomial trends (3 per component)
[..+nS-1]                             | Background RBF weights
[..+3-1]                              | Background polynomial trend (3 coeffs)
[..+2*nS]                             | Stamp positions (x₀, y₀, x₁, y₁, ...)
```

**Global State Parameters:**

- `useTPS` (int, default D_USE_TPS): Enable TPS mode (1) or polynomial (0)
- `tpsSmoothing` (double, default D_TPS_SMOOTHING): Regularization parameter (1e-6)
- `nS` (int): Number of stamps per region (used to calculate offsets)

**Python API Integration:**

```python
config = KernelConfig(
    use_tps=True,              # Enable TPS mode
    tps_smoothing=1e-6         # Regularization strength
)
kernel_sol = fit_kernel(science, template, config=config)
output = spatial_convolve(template, kernel_sol, config=config)
```

**CLI Usage:**

```bash
./hotpants -iimage.fits -ttemplate.fits -ooutput.fits \
  -nrx 1 -nry 1 -nsx 10 -nsy 10 \
  -useTPS 1 -tpsSmoothing 1e-6
```

**Performance Characteristics:**

- **Fitting overhead:** ~2–5% additional time (RBF matrix assembly + LAPACK LU solve)
- **Evaluation overhead:** Negligible (~1 FLOP per kernel + background evaluation)
- **Memory overhead:** ~(nS + 5) doubles per component (for RBF weights + poly trend + positions)
- **Benefit:** Enables single-region full-frame processing (no tiling artifacts)

**Error Handling:**

- If `tps_fit_kernel()` fails (singular matrix), falls back to polynomial; logs error
- If `tps_fit_background()` fails, background reverts to polynomial; logs warning
- Graceful degradation ensures robustness on difficult datasets

**Future Extensions:**

- Regularization options (L2, Laplacian, etc.) beyond current single parameter
- Adaptive smoothing based on data quality or residual statistics
- Multi-resolution TPS (coarse grid with refinement) for very large images
- Integration with Bramich (2008) delta-function basis (see Delta Basis section below)

---

### Delta Function Kernel Basis Implementation (Phases 1–6, May 2026)

**Status:** Phases 1–6 complete; pluggable infrastructure and core algorithms working.
Phases 3+ and 5+ (extended work) remain for full fitKernel integration.

**Overview:**

The Bramich (2008) delta function basis provides an alternative to the multi-Gaussian
basis for kernel representation. Instead of a fixed set of Gaussian profiles weighted
by spatial polynomials, delta basis uses an adaptive grid of narrow RBF functions
(soft deltas via thin-plate spline kernel φ(r) = r² log(r)), offering:

- More direct kernel representation (directly solve for grid values)
- Flexible spatial resolution (via `deltaKerGridSize` parameter)
- Automatic smoothness control (Laplacian regularization)
- Foundation for future basis variants (PCA, Fourier, etc.)

**Architecture:**

Delta basis uses a pluggable architecture with clean dispatch points:

1. **Basis Type Enumeration:** `iBasisType` global (0=Gaussian, 1=Delta)
   - Routed at: `getKernelVec()`, `make_kernel_dispatch()`, `spatial_convolve()`
   - Existing Gaussian path unchanged; delta code isolated

2. **Grid Context:** `delta_basis_grid_t` struct stores grid geometry
   - `ngrid_x, ngrid_y`: Grid dimensions
   - `nbasis`: Total basis functions (product of grid dims)
   - `grid_cx[], grid_cy[]`: RBF center coordinates
   - `grid_spacing`: User-specified distance between points (pixels)

3. **RBF Evaluation:**
   - `delta_rbf_kernel(r)`: Thin-plate spline RBF φ(r) = r² log(r)
   - `eval_delta_basis(idx, dx, dy)`: Evaluate basis function i at offset (dx, dy)
   - Numerically stable (r < ZEROVAL → φ(r) = 0)

4. **Stamp Convolution:**
   - `xy_conv_stamp_delta()`: Direct 2D correlation of stamp with delta RBF
   - O(k²m²) complexity (k=kernel width, m=stamp width)
   - Evaluates RBF on-the-fly (no pre-computed arrays)
   - Zero-padded boundary handling

5. **Regularization:**
   - `assemble_laplacian_regularization()`: Builds discrete 2D Laplacian matrix
   - 5-point stencil (cardinal neighbors only)
   - Computes R = λ·Lᵀ·L added to normal equations
   - Prevents overfitting and noise amplification

6. **Kernel Evaluation & Convolution:**
   - `make_kernel_delta()`: Evaluate spatially-varying kernel [STUB in Phase 5]
   - `spatial_convolve_delta()`: FFT-accelerated convolution [STUB in Phase 5]
   - Integrates with polynomial/TPS spatial variation

**Configuration Parameters:**

| Global | Type | CLI Flag | Python | Default | Range |
|--------|------|----------|--------|---------|-------|
| `iBasisType` | int | `-basisType` | `basis_type` | 0 (gaussian) | 0, 1 |
| `rDeltaKerGridSize` | double | `-deltaKerGridSize` | `delta_ker_grid_size` | 2.0 px | > 0 |
| `rDeltaRegularization` | double | `-deltaRegularization` | `delta_regularization` | 1e-3 | ≥ 0 |

**Grid Initialization (`init_delta_basis_grid`):**

1. Validate hwKernel and grid spacing
2. Compute grid dimensions: ngrid_x = ceil(2·hwKernel / spacing) + 1
3. Allocate grid_cx[], grid_cy[]
4. Populate grid points from -hwKernel to +hwKernel at specified spacing
5. Update global nCompKer = nbasis

**Example Grid (hwKernel=10, spacing=2.0):**

```
Grid dims: 11×11 = 121 basis functions
X: [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
Y: [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
```

**Regularization Strategy:**

Laplacian penalty minimizes kernel curvature:

```
min ||A·x - b||² + λ·||L·x||²
```

where L is the discrete 2D Laplacian. Effect:

- λ = 0: No regularization (noisy kernel, worse fit stability)
- λ = 1e-3: Mild smoothing (good balance)
- λ > 1e-2: Strong smoothing (may degrade fit quality)

**Integration Points:**

- `getKernelVec()`: Calls `init_delta_basis_grid()` when iBasisType=1
- `fillStamp()`: Routes to `xy_conv_stamp_delta()` per basis type [Phase 3+]
- `fitKernel()`: Applies `assemble_laplacian_regularization()` before LAPACK solve
- `make_kernel_dispatch()`: Routes to `make_kernel_delta()` per basis type
- `spatial_convolve()`: Routes to `spatial_convolve_delta()` per basis type

**Python API Usage:**

```python
from hotpants import fit_kernel, spatial_convolve, KernelConfig

# Configure delta basis
config = KernelConfig(
    basis_type="delta",              # Use delta RBF basis
    delta_ker_grid_size=2.0,        # 2 pixel spacing
    delta_regularization=1e-3,      # Smoothness penalty
    use_tps=False                   # Or True for TPS spatial variation
)

# Fit kernel
kernel_sol = fit_kernel(template, science, config=config)

# Apply convolution
diff = spatial_convolve(science, kernel_sol, config=config)
```

**Current Limitations (Phase 7 status):**

1. `build_matrix_delta()` not fully implemented (normal equation accumulation)
   - Requires integration with fillStamp() and stamp structure
   - Deferred to extended work (Phase 3+)

2. `make_kernel_delta()` and `spatial_convolve_delta()` are stubs
   - Interface defined; full implementation deferred (Phase 5+)
   - Requires careful integration with spatial variation logic

3. No unit tests yet for delta basis (framework ready, tests pending)

4. No performance optimization (code correct but may have marginal inefficiencies)

**Future Work:**

- Complete `build_matrix_delta()` and integrate with fitKernel() pipeline
- Implement `make_kernel_delta()` with spatial variation support
- Implement `spatial_convolve_delta()` with FFT acceleration
- Add comprehensive unit and integration tests
- Performance profiling and optimization
- Support for additional basis types (PCA, Fourier, etc.) via same pluggable framework

**Reference:** Bramich (2008), "The Optimal Difference Image Combination of Dithered Images",
MNRAS 389:1365 (arXiv:0802.1273)

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

### ✓ COMPLETED: Thin Plate Splines for Spatial Variation

**Status:** Fully implemented and tested (May 2026)

**Implementation Details:**
- Thin Plate Spline (TPS) RBF surfaces replace polynomial spatial variation for kernel
  coefficients and background
- Enables smooth full-frame processing without tiling artifacts (use `-nrx 1 -nry 1`)
- Accessible via Python API: `KernelConfig(use_tps=True, tps_smoothing=1e-6)`
- CLI flags: `-useTPS 1 -tpsSmoothing <value>`
- 84/84 tests passing; graceful fallback to polynomial if TPS fitting fails

**Key Functions Added:**
- `tps_kernel()`, `tps_assemble_matrix()`, `tps_fit_coefficients()`, `tps_evaluate()`
- `tps_fit_kernel()`, `make_kernel_tps()`, `tps_fit_background()`, `get_background_tps()`

**Benefits:**
- Smooth spatial variation without polynomial step discontinuities
- Enables practical single-region full-frame processing for widefield instruments
- Regularization parameter (`tps_smoothing`) controls overfitting vs. fit quality

---

* ✓ **PARTIAL: Bramich (2008) Delta Function Kernel Basis** (May 2026, Phases 1-6 complete)
  
  **Status:** Pluggable basis infrastructure complete. Core grid initialization, RBF
  evaluation, and stamp convolution working. Laplacian regularization implemented.
  Kernel evaluation and full convolution stubs in place. Python API fully integrated.
  
  **What's done:**
  - Phase 1: Basis type enum, CLI flags, C API wrapper functions
  - Phase 2: Delta grid initialization with configurable RBF spacing
  - Phase 3: Direct 2D convolution with delta RBF basis (xy_conv_stamp_delta)
  - Phase 4: Discrete 2D Laplacian regularization matrix assembly
  - Phase 5: Kernel evaluation and convolution function stubs (interface defined)
  - Phase 6: Full Python API support (KernelConfig.basis_type, delta parameters)
  
  **What remains:**
  - Phase 3+ (extended): Full normal equation accumulation (build_matrix_delta)
  - Phase 5+ (extended): Complete kernel evaluation with spatial variation
  - Phase 5+ (extended): Full FFT-based convolution with delta kernels
  - Integration with fitKernel() pipeline and spatial variation (polynomial/TPS)
  
  **Usage (Python API):**
  ```python
  config = KernelConfig(
      basis_type="delta",           # Use delta RBF basis
      delta_ker_grid_size=2.0,     # Grid spacing (pixels)
      delta_regularization=1e-3    # Laplacian smoothness penalty
  )
  solution = fit_kernel(template, science, config=config)
  ```
  
  **CLI Usage:**
  ```bash
  hotpants -basisType 1 -deltaKerGridSize 2.0 -deltaRegularization 1e-3 ...
  ```
  
  **Architecture:** Pluggable basis system allows future extensions (PCA, Fourier, etc.)
  without modifying core Gaussian code path. All dispatch points implemented with
  clear separation of concerns.

### ✓ COMPLETED: LSST DMTN-021 Afterburner Decorrelation (May 2026)

**Status:** Fully implemented with spatially-varying whitening kernel and FITS output

**Algorithm Overview:**

The DMTN-021 decorrelation removes correlated noise introduced by PSF-matching convolution
using a spatially-varying whitening kernel φ(x,y):

```
φ̂(k) = √[(σ_s² + σ_t²) / (σ_s² + σ_t² · κ̂²(k))]
D_decorr = φ ⊗ D
```

where κ̂(k) is the fitted kernel in Fourier space, σ_s² and σ_t² are science and template
image variances, and ⊗ denotes spatial convolution.

**Implementation Details:**

1. **Core Decorrelation (decorrelation.c, ~750 lines):**
   - `compute_decorrelation_kernel_fft()`: Applies DMTN-021 formula with numerical stability
   - `fit_decorrelation_grid()`: Evaluates φ at stamp centers via local kernel reconstruction
   - `eval_decorrelation_at_point()`: Bilinear interpolation of φ field to image pixels
   - `spatial_convolve_with_decorr()`: Applies φ ⊗ D via direct spatial convolution
   - `decorrelation_init_region()`: High-level wrapper after kernel fitting
   - `decorrelation_apply_region()`: High-level wrapper with float/double marshalling

2. **Integration Points (main.c):**
   - Decorrelation triggered after `fitKernel()` per region (line ~1206)
   - Decorrelated image computed via `decorrelation_apply_region()` (line ~1280)
   - Output written to FITS layer 2 with metadata keywords (line ~1768)
   - FITS EXTNAME='DECORR' with header keywords for method and variances

3. **CLI Configuration (vargs.c):**
   - `-useDecorrelation [0|1]`: Enable/disable decorrelation (default: 0)
   - `-decorrScienceVar <float>`: Science image variance (default: auto-estimate)
   - `-decorrTemplateVar <float>`: Template image variance (default: auto-estimate)
   - `-decorrUseTPS [0|1]`: Use TPS interpolation for φ field (default: 0, bilinear)

4. **Numerical Stability:**
   - Avoids division-by-zero when denominator is very small (φ_min = 1e-6)
   - Clamps φ to physically meaningful range [1e-6, 1e6]
   - Handles zero-kernel case gracefully (φ >> 1, amplifies noise)
   - Works correctly with small and large variances (tested up to 1e-10 and 1e10)

5. **Test Coverage (test_decorrelation.py, 19 tests, all passing):**
   - Configuration and formula validation
   - Bilinear interpolation correctness and monotonicity
   - Decorrelation effects on synthetic data (uniform and spatially-varying)
   - Integration with kernel fitting pipeline
   - Numerical stability with extreme variances
   - Edge cases (1×1 images, single row/column, zero images)

**Usage:**

CLI:
```bash
./build/hotpants -iimage.fits -ttemplate.fits -ooutput.fits \
  -useDecorrelation 1 -decorrScienceVar 100.0 -decorrTemplateVar 100.0
```

Output FITS structure:
```
Primary HDU
├─ Data: Original difference image D
├─ Keywords: EXTNAME='DIFF'
└─ Comments: Decorrelation method, variances used

HDU 2: Decorrelated difference image D_decorr = φ ⊗ D
├─ Data: D_decorr
├─ EXTNAME: 'DECORR'
├─ HIERARCH DECORR_METHOD: 'AFTERBURNER'
├─ HIERARCH DECORR_SCI_VAR: 100.0
└─ HIERARCH DECORR_TPL_VAR: 100.0
```

**Variance Estimation:**

If not provided via CLI flags, variances are auto-estimated from image statistics:
- Science variance: σ-clipped statistics of science image
- Template variance: σ-clipped statistics of template image

**Performance:**

- Decorrelation grid fitting: ~1–5% additional time per region (FFT of kernel + grid eval)
- Decorrelated convolution: ~1–2× FFT convolution cost (spatially-varying φ requires per-pixel eval)
- Memory overhead: O(nStamps) for grid storage (typically <1 MB)

**Known Limitations:**

1. Bilinear interpolation is the default; TPS interpolation planned as future optimization
2. Global variances used; future work could support per-region variance estimation
3. Decorrelation disabled in multi-region mode by default (designed for single-region with TPS)

**Future Extensions:**

- TPS interpolation of φ field for smoother spatial variation
- Per-region variance estimation for heterogeneous image statistics
- Integration with noise covariance estimation from residual images
- Support for cascaded decorrelation (iterative refinement)

**References:**

- LSST Data Management Technical Note DMTN-021: "Image Difference Decorrelation"
  https://dmtn-021.lsst.io/
- Details on formula, numerical stability, and spatially-varying implementation

---

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
