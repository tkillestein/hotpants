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
├── main.c              # Pipeline orchestration: FITS I/O, region loop, output assembly
├── alard.c             # Core algorithm: kernel fitting, spatial convolution, LU solver
├── functions.c         # Stamps, PSF-centre finding, statistics, masking
├── vargs.c             # CLI argument parsing (~50 options)
├── defaults.h          # Compile-time parameter defaults
├── globals.h           # Global variable declarations
├── functions.h         # Function prototypes and CFITSIO includes
├── extractkern.c       # Standalone utility: reconstruct and inspect kernel
├── extractkernOnes.c   # Variant of extractkern
├── maskim.c            # Standalone utility: apply mask to image
├── Makefile            # Linux build
├── Makefile.macosx     # macOS build
├── README.md           # User-facing option reference
├── NOTES               # Version history and profiling data
└── PROF                # gprof output from representative runs
```

### Key functions

| Function | File | Purpose |
|---|---|---|
| `main()` | `main.c` | Orchestrates the full pipeline |
| `fitKernel()` | `alard.c` | Least-squares kernel fit per region |
| `build_matrix()` / `build_scprod()` | `alard.c` | Accumulate normal equations |
| `ludcmp()` / `lubksb()` | `alard.c` | LU decomposition and back-substitution (Numerical Recipes) |
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
# Linux
make hotpants

# macOS
make -f Makefile.macosx hotpants
```

**Compiler flags:** `-O3 -funroll-loops -std=c99 -pedantic-errors -Wall`

**External dependency:** CFITSIO (`-lcfitsio`). Edit `CFITSIOINCDIR` and
`LIBDIR` in the Makefile to match your installation.

No test suite or CI exists yet.

---

## Performance Profile

From `PROF` and `NOTES` (gprof on representative runs):

| Function | Approx. CPU share |
|---|---|
| `spatial_convolve()` | ~60% |
| `getPsfCenters()` | ~35% |
| `xy_conv_stamp()` | ~18% (subset of above) |
| LU decomposition | < 5% |

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

### Parallelism — OpenMP and FFT convolution

**OpenMP**

- The region loop in `main.c` is the natural top-level parallelism boundary:
  regions are fully independent after the initial FITS read.
- Add `#pragma omp parallel for schedule(dynamic)` to the region loop.
- `spatial_convolve()` (`alard.c`) is embarrassingly parallel at the pixel
  level; the `j1` outer loop can be parallelised with OpenMP.
- Shared mutable state: `mRData` mask array is updated during PSF-centre
  finding — protect with `#pragma omp critical` or restructure to local buffers
  merged after the parallel section.
- Compile with `-fopenmp`; respect `OMP_NUM_THREADS` at runtime.
- The previous thread-pool approach (`thpool.c`/`thpool.h`) has been removed
  from this branch; do not reintroduce it — OpenMP is sufficient and simpler.

**FFT convolution (high-value optimisation)**

The kernel decomposition K = Σᵢ cᵢ(x,y)·φᵢ means the spatially-varying
convolution can be restructured as:

```
D(x,y) = Σᵢ cᵢ(x,y) · [I ⊗ φᵢ](x,y)
```

Each I⊗φᵢ is a convolution with a **fixed** kernel, computable in O(N log N)
via FFT. The final image is then a pixel-wise weighted sum of N_basis ≈ 50
precomputed images. This replaces the O(N·k²) direct loop in
`spatial_convolve()` and should give substantial speedups for large images.

Implement via **FFTW3** (`-lfftw3`). Plan the FFTs once per region and reuse
plans across basis elements.

### Build system

- Migrate from `Makefile` to **CMake** to cleanly handle the new optional
  dependencies (LAPACK/BLAS, FFTW3, OpenMP) and the Python extension build.
- Keep a thin `Makefile` shim at the root for convenience (`make` → `cmake
  --build`).
