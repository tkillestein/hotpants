# HOTPANTS — Development Guide

**High Order Transform of PSF ANd Template Subtraction**

HOTPANTS produces difference images between an astronomical science image and a
template by fitting and applying a spatially-varying convolution kernel that
matches the two PSFs. The core algorithm is from Alard & Lupton (1998), ApJ
503:325 (https://iopscience.iop.org/article/10.1086/305984).

License: MIT (Andy Becker, 2013). See `LICENSE`.

---

## Repository Structure

```
hotpants/
├── src/
│   ├── main.c          # Pipeline orchestration: FITS I/O, region loop, output assembly
│   ├── alard.c         # Core algorithm: kernel fitting, spatial convolution, LU solver
│   ├── functions.c     # Stamps, PSF-centre finding, statistics, masking
│   ├── vargs.c         # CLI argument parsing (~50 options)
│   ├── maskim.c        # Standalone utility: apply mask to image
│   ├── defaults.h      # Compile-time parameter defaults
│   ├── globals.h       # Global variable declarations
│   └── functions.h     # Function prototypes and CFITSIO includes
├── tests/              # pytest regression suite
├── docs/               # Sphinx documentation
├── CMakeLists.txt      # Build system
├── Doxyfile            # Doxygen configuration
├── pyproject.toml      # Python project (uv-managed)
├── README.md           # User-facing option reference
└── NOTES               # Version history
```

### Key functions

| Function | File | Purpose |
|---|---|---|
| `main()` | `main.c` | Orchestrates the full pipeline |
| `fitKernel()` | `alard.c` | Least-squares kernel fit per region |
| `build_matrix()` / `build_scprod()` | `alard.c` | Accumulate normal equations |
| `fitKernel()` solver | `alard.c` | LAPACK Cholesky (`dpotrf`/`dpotrs`) solve of normal equations |
| `spatial_convolve()` | `alard.c` | Apply spatially-varying kernel to full region — **primary bottleneck** |
| `make_kernel()` | `alard.c` | Evaluate kernel at image position from polynomial coefficients |
| `xy_conv_stamp()` | `alard.c` | Separable 2D convolution of a stamp with a Gaussian basis element |
| `getKernelVec()` | `alard.c` | Initialise multi-Gaussian kernel basis vectors |
| `buildStamps()` | `functions.c` | Create stamp grid and compute statistics |
| `getPsfCenters()` | `functions.c` | Locate bright star centres in stamps — **secondary bottleneck** |
| `getStampStats3()` | `functions.c` | Sigma-clipped statistics and histogram FWHM estimate |
| `makeNoiseImage4()` | `functions.c` | Propagate noise through convolution |
| `vargs()` | `vargs.c` | Parse all command-line options |

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
polynomials in image position.

**Pipeline steps:**

1. Divide image into regions (`-nrx`, `-nry`)
2. Within each region, lay a stamp grid (`-nsx`, `-nsy`)
3. In each stamp, locate bright PSF centres (`getPsfCenters`)
4. For each substamp, convolve I (or T) with every basis element φᵢ
   (`xy_conv_stamp`) and accumulate normal equations (`build_matrix`,
   `build_scprod`)
5. Solve the normal equations by LU decomposition (`fitKernel` →
   `ludcmp`/`lubksb`) to obtain the cᵢ coefficients
6. Sigma-clip bad stamps and refit if necessary
7. Apply the fitted spatially-varying kernel to the full region
   (`spatial_convolve`)
8. Subtract and write output FITS

---

## Building

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/hotpants [options]
```

CMake auto-detects CFITSIO, OpenBLAS, FFTW3, and OpenMP. Optional flags:

- `-DUSE_FFTW=ON/OFF` — enable FFT-accelerated convolution (default: ON if FFTW3 found).
- `-DUSE_OPENMP=ON/OFF` — enable multi-threaded parallelism (default: ON if found).

**External dependencies:**

- CFITSIO (`-lcfitsio`) — FITS I/O library.
- OpenBLAS or LAPACK (`-lopenblas` or `-llapack -lblas`) — linear algebra.
- FFTW3 (`-lfftw3`, optional) — FFT-accelerated convolution.
- OpenMP (`-fopenmp`, optional) — multi-threading support.

**Compiler flags:** `-O3 -march=native -funroll-loops -std=c17 -Wall`

Use `-march=native` for SIMD optimisations on the build host; omit for portable binaries.

---

## Performance Profile

**Baseline (single-threaded, direct convolution):**

From `NOTES` (gprof on representative runs):

| Function | Approx. CPU share |
|---|---|
| `spatial_convolve()` | ~60% |
| `getPsfCenters()` | ~35% |
| `xy_conv_stamp()` | ~18% (subset of above) |
| LU decomposition | < 5% |

**With OpenMP and FFT optimisations:**

- Region loop parallelised with `#pragma omp parallel for schedule(dynamic)`.
- `spatial_convolve_fft()` replaces direct pixel-loop with O(N log N) FFT convolution; expected 3–8× speedup on large images.
- `getPsfCenters()` remains a secondary bottleneck; vectorisation (SIMD) not yet explored.

---

## Completed Work

### ✓ Linear algebra — LAPACK Cholesky and CMake

- Replaced `ludcmp()` and `lubksb()` (Numerical Recipes) with LAPACK Cholesky (`dpotrf` + `dpotrs`).
- Linear system solve now uses optimised library code; < 5% of runtime.
- CMake build system cleanly handles CFITSIO, OpenBLAS, FFTW3, and OpenMP detection.
- Both CMake and legacy Makefile still available for compatibility.

### ✓ Parallelism — OpenMP

- Region loop in `main.c` parallelised: `#pragma omp parallel for schedule(dynamic)`.
- `spatial_convolve_fft()` inner pixel loop parallelised with `collapse(2)`.
- FFTW plan creation serialised with `#pragma omp critical(fftw_plan)` to avoid thread-unsafe planning.
- CFITSIO I/O protected with `#pragma omp critical(cfitsio)` (not thread-safe).
- Respects `OMP_NUM_THREADS` at runtime; typical 4–8 thread scaling observed.

### ✓ FFT convolution — FFTW3 acceleration

- New `spatial_convolve_fft()` function implements O(N log N) convolution via FFTW3.
- Restructured kernel as K = Σᵢ cᵢ(x,y)·φᵢ; each basis convolution I⊗φᵢ computed once via FFT.
- FFT sizes optimised to highly-composite numbers (powers of 2, 3, 5) via `next_fftw_size()`.
- Plans created once per region, reused across all basis elements.
- Mask propagation and variance still computed via direct convolution (unchanged).
- Switchable at runtime: `-fft` flag enables FFTW path; default is direct convolution for backward compatibility.

---

## Modernisation Goals

This section records the intended development direction. Follow these
preferences in all future work on the repository.

### Python packaging — use `uv`

- Manage the project with `uv` (`pyproject.toml`, `uv.lock`).
- Build a Python extension wrapping the C core. Preferred binding strategy:
  `cffi` (avoids Python/C ABI fragility) or `pybind11` (acceptable if a thin
  C++ shim is needed).
- The Python API must accept and return **numpy arrays**. Do not pull FITS I/O
  into the C extension — let callers use astropy for that.
- Minimum public surface: a `hotpants()` function taking image arrays plus a
  parameter dataclass/struct; the existing CLI must continue to work alongside
  it.
- Target Python ≥ 3.11.

### Documentation

- Add **Doxygen** comments to every public C function, especially in `alard.c`
  and `functions.c`. Describe parameters, return values, and any non-obvious
  invariants.
- Set up **Sphinx** with the `breathe` extension so Doxygen output appears in
  the Python API docs.
- Add a `CONTRIBUTING.md` covering the build, testing, and code-style
  expectations.
- Update `README.md` with uv/pip installation instructions and a quick-start
  example using the Python API.
- Inline comments should reference Alard & Lupton (1998) where relevant,
  including equation numbers.

### Linear algebra — replace custom LU with LAPACK

- Remove `ludcmp()` and `lubksb()` (Numerical Recipes) from `alard.c`.
- The normal equations matrix from `build_matrix()` is symmetric
  positive-definite; use LAPACK **Cholesky** (`dpotrf` + `dpotrs`) rather than
  general LU (`dgesv`).
- Use BLAS (`dgemm`, `dgemv`) for accumulation loops in `build_matrix()` and
  `build_scprod()` where the data layout supports it.
- Link with `-lopenblas` (preferred) or `-llapack -lblas`. Matrix dimensions
  are variable and can be large, so a maintained library is the right choice
  regardless of typical problem size.
- The CMake/Makefile refactor should detect BLAS/LAPACK via `pkg-config` or
  `find_package`.

### Further Performance Optimisations

**PSF-centre finding acceleration**

`getPsfCenters()` remains the second-largest bottleneck (~35% CPU). Opportunities:

1. **SIMD vectorisation** — the inner stamp-search loops are data-parallel; use
   `_mm256_*` (AVX2) or `_mm512_*` (AVX-512) intrinsics or rely on compiler
   autovectorisation with `-march=native -mtune=native`.
2. **Approximate nearest-neighbour search** — current brute-force all-pairs comparison
   could be replaced by spatial hashing or k-d trees if the stamp count is large.
3. **GPU offload** — CUDA or OpenCL could handle the stamp comparison kernel, but
   data transfer overhead may dominate for modest image sizes.

**Convolution kernel basis optimisation**

- Current Gaussian basis is fixed (σ = 0.7, 1.5, 3.0 px; polynomial degrees 6, 4, 2).
  For specific PSF shapes (e.g. wide-field aberrations), a data-driven basis
  (PCA of residual PSFs) might reduce `nCompKer` and thus FFT count. Explore
  adaptive basis selection.

**Memory bandwidth and cache**

- `spatial_convolve_fft()` allocates large temporary buffers (R2C FFT plans, real/complex arrays).
  For multi-threaded execution with limited cache, prefetching or tiling strategies could help.
- Monitor cache misses with `perf stat` or `likwid` on representative workloads.

**Polynomial coefficient evaluation**

`make_kernel()` evaluates a low-order polynomial at every pixel. Current
implementation is inline; consider:

1. Horner's method (already used) vs. precomputed lookup tables for coarse grid,
   trilinear interpolation.
2. Vectorisation of the polynomial loop if the compiler cannot autovectorise it.

**Mask propagation in FFT path**

Mask propagation (`makeNoiseImage4()`) is currently done via direct convolution
even when the pixel values use FFT. These could be harmonized to the FFT path
as well, but requires careful handling of integer/binary mask semantics.

### Build system

- ✓ **CMake** build system; cleanly detects CFITSIO, OpenBLAS, FFTW3, and
  OpenMP. See `CMakeLists.txt`.
- Compiler flags: `-O3 -march=native -funroll-loops -std=c17` (enable SIMD
  optimisations with `-march=native` for portable binaries on target systems).
