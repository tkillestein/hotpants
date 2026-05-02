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

---

## Readability Improvements

This section identifies code readability antipatterns discovered in the HOTPANTS
codebase and proposes improvements to make the implementation clearer, more
maintainable, and faster to onboard new contributors (especially those building
the Python API). Clear code is essential for debugging, optimization, and
understanding the complex algorithms from Alard & Lupton (1998).

### Readability Antipatterns

The table below catalogues specific readability issues, their locations,
impact, and recommended fixes. Use this as a reference guide when refactoring
or writing new code.

| Antipattern | Category | Found in | Problem | Fix | Priority | Effort |
|---|---|---|---|---|---|---|
| Global variable prefix abbreviations (tRData, iRData, mRData, oRData, eRData) | Naming | globals.h:37-75 | No prefix legend; forces reader to deduce meaning from context | Add comment block in globals.h explaining: `t*` = template, `i*` = image, `m*` = mask, `o*` = output, `e*` = ephemeris | HIGH | 30 min |
| Single-letter loop variables in nested contexts (i, j, k, l, m, n) | Naming | alard.c, functions.c (pervasive) | Confusing when nesting depth > 2; hard to trace what variable refers to | Use descriptive names: `pixel_x`, `pixel_y`, `bin_idx`, `vector_idx`, `stamp_idx` based on loop purpose | HIGH | 2-3 hrs |
| Cryptic flag names: `ren` | Naming | alard.c:108, 111, 116 | Purpose obscured; requires tracing code to understand as "renormalization flag" | Rename to `renormalize_flag` and document when/why it's set | MEDIUM | 1 hr |
| Cryptic variable names: `q`, `ncomp1`, `ncomp2`, `ivecbg` | Naming | alard.c:589, 577-578, 571 | Intent hidden; `q` = scalar products in normal equations, others relate to polynomial basis structure | Use full names: `scalar_product`, `num_kernel_components_low`, `num_kernel_components_high`, `background_vector_index` | MEDIUM | 1-2 hrs |
| Mixed naming conventions (camelCase vs snake_case) | Naming | globals.h, alard.c (inconsistent) | Legacy code uses `kerOrder`, `fwKernel`; newer uses `build_matrix`, `spatial_convolve` | Document preference (snake_case for new code, accept camelCase for legacy globals); add convention rule to contributor checklist | LOW | Ongoing |
| Magic formula: (kerOrder+1)*(kerOrder+2)/2 | Magic Numbers | alard.c:499, 577, 640, 726, etc. (9+ occurrences) | Formula appears repeatedly without explanation; connection to triangular polynomial basis not documented | Replace with named constant `NUM_POLY_TERMS(order)` or `#define` and add comment: "Triangular basis: polynomial degrees i,j ≤ order with i+j ≤ order; term count = (order+1)*(order+2)/2. See Alard & Lupton 1998, Eq. 3." | HIGH | 1-2 hrs |
| Parity check arithmetic: (idegx / 2) * 2 - idegx | Magic Numbers | alard.c:109-110 | Computes odd/even detection without comment; cryptic bit-twiddling pattern | Add comment explaining: "Detects odd-degree terms; result is 1 if idegx is odd, 0 if even. Used to trigger renormalization for higher-order basis functions." | MEDIUM | 30 min |
| Hardcoded histogram thresholds: nstat=100, ufstat=0.9, mfstat=0.5 | Magic Numbers | functions.c:1008-1010 | Why these specific values? Rationale for 100 pixels, 0.9 and 0.5 percentiles not explained | Define in defaults.h with names `HISTOGRAM_SAMPLE_SIZE`, `HISTOGRAM_UPPER_FRAC`, `HISTOGRAM_LOWER_FRAC` and add comments explaining histogram bin-width estimation | HIGH | 1 hr |
| Hardcoded RNG seed: idum = -666 | Magic Numbers | functions.c:1023 | "The devil's seed" comment mentions Numerical Recipes convention but algorithm rationale absent | Define `RNG_SEED_MAGIC -666` in defaults.h with comment: "Numerical Recipes convention; triggers internal re-seeding of Park & Miller generator" | MEDIUM | 30 min |
| Sigma-clipping iteration cap: tries >= 5 with no explanation | Magic Numbers | functions.c:1094 | Why 5? What does iteration do? Convergence criteria not documented | Define `MAX_SIGMA_CLIP_RETRIES 5` in defaults.h with comment: "Iteration cap to prevent divergence in histogram-based bin-width estimation. Typical convergence < 3 iterations." | MEDIUM | 30 min |
| Undocumented polynomial basis construction | Documentation | alard.c:104-120 | Triple nested loop (ig → idegx → idegy) with variable `nvec` counter; relationship between triplet and vector index is implicit | Add Doxygen comment before loop explaining: "Triangular basis construction: for each Gaussian basis element φᵢ (ig), iterate over polynomial degrees (idegx, idegy) such that idegx+idegy ≤ kerOrder. Assemble convolution responses for each (ig, idegx, idegy) triplet into vector nvec. See Alard & Lupton 1998, Sect. 2.2." | HIGH | 1 hr |
| Undocumented separable convolution intent | Documentation | alard.c:294-358 | Comment mentions "separability" but doesn't explain why it's chosen or the performance benefit | Add comment: "Gaussian-polynomial filters are separable: G(x,y)*P(x,y) = [G(x)*P_x(x)] * [G(y)*P_y(y)]. This reduces 2D convolution O(k²n²) to two 1D passes O(kn²), critical for performance. Renormalization flag adjusts for half-Gaussian overlap." | MEDIUM | 30 min |
| Undocumented histogram FWHM algorithm | Documentation | functions.c:979-1100+ | getStampStats3() performs 5-iteration loop for histogram bin refinement; nested loops, thresholding logic unmarked | Extract algorithm into separate function with Doxygen block: "Algorithm: (1) sample 100 random pixels to estimate range, (2) compute bin width from percentile ranges, (3) build 256-bin histogram, (4) locate peak (mode) via threshold, (5) iterate if necessary for convergence. Returns FWHM and mean." | MEDIUM | 2 hrs |
| Oversized functions without segmentation (check_stamps, spatial_convolve_fft, getStampStats3) | Code Structure | alard.c:705-960, alard.c:1684+, functions.c:979+ | Functions exceed 100-150 lines; multiple concerns intermingled (matrix extraction, solving, statistics); hard to test/debug individual steps | Break into smaller sub-functions: `extract_stamp_matrix()`, `solve_kernel_linear_system()`, `compute_sigma_clipped_mean()`, `estimate_histogram_fwhm()`. Document with Doxygen. | MEDIUM | 4-6 hrs |
| Scattered obsolete/commented-out debug code | Code Quality | functions.c:270, 1095 | Fragments like `/* DD fprintf(...) */` and commented variable definitions create noise; unclear why kept | Remove dead code entirely; use version control for history. Retain only active debug output with explicit `verbose` level checks. | LOW | 1-2 hrs |
| Inconsistent error handling and logging | Code Quality | alard.c:92-100, functions.c:1095-1098 | Mix of `verbose >= 1` and `verbose >= 2` checks; ad-hoc `fprintf` calls without standardized format | Implement centralized logging macro: `VERBOSE_LOG(level, fmt, ...)` that checks verbosity level and prints to stderr with consistent formatting. Document verbose levels (0=silent, 1=progress, 2=debug). | LOW | 2-3 hrs |

### Improvement Priorities & Roadmap

**HIGH PRIORITY (Quick Wins, Direct Impact):** 1–3 hours per item

- Add global variable naming legend to globals.h (30 min)
- Replace magic numbers with #define constants in defaults.h (1-2 hrs)
- Document polynomial basis formula in comments (1 hr)
- Add Doxygen summaries to complex functions (2-3 hrs)

**MEDIUM PRIORITY (Refactoring, Medium Effort):** 1–6 hours per item

- Rename loop variables in complex nesting contexts (2-3 hrs per file)
- Extract sub-functions from getStampStats3 to clarify histogram algorithm (2 hrs)
- Document parity checks and renormalization logic (1 hr)
- Consolidate logging/verbosity checks (2-3 hrs)

**LOW PRIORITY (Nice-to-Have, Larger Refactoring):** 4–10+ hours per item

- Full function decomposition of spatial_convolve_fft and xy_conv_stamp (4-6 hrs each)
- Standardize naming conventions across all files (8-10 hrs)
- Remove all dead code and obsolete comments (1-2 hrs)

### Actionable Guidance by Category

**A. Global Variable Naming**

Add a naming legend at the top of `globals.h`, right after the struct definitions:

```
/*
  GLOBAL VARIABLE PREFIX LEGEND
  ==============================
  These variables are declared with single-letter prefixes for brevity but
  clarity is essential. Prefix meanings:
  
  t* = template image (e.g., tGain, tUThresh, tRData)
  i* = science image (e.g., iGain, iUKThresh, iRData)  
  m* = mask/masking data (e.g., mRData, mGain)
  o* = output image (e.g., outim, oRData)
  e* = reserved for ephemeris/environment (currently unused)
  
  Example: tUThresh = template upper threshold for saturation
           iRdnoise = science image readnoise
*/
```

**B. Loop Variables in Complex Contexts**

When nesting depth exceeds 2 or when loop purpose isn't immediately obvious from context,
use descriptive names. Apply this rule to files like `alard.c` (kernel basis iteration,
polynomial assembly) and `functions.c` (histogram binning, stamp processing).

Rule of thumb: If you cannot immediately name the loop variable's domain (e.g., "pixel
columns", "polynomial degrees", "histogram bins"), it's too cryptic and needs renaming.

**C. Magic Numbers Conversion**

Convert all "magic" constants to named `#define` in `defaults.h` with explanatory comments:

- `HISTOGRAM_SAMPLE_SIZE 100` — number of random pixels sampled to estimate bin width
- `HISTOGRAM_NUM_BINS 256` — resolution of final histogram for FWHM estimation
- `HISTOGRAM_LOWER_FRAC 0.5` — lower percentile for bin-width estimation
- `HISTOGRAM_UPPER_FRAC 0.9` — upper percentile for bin-width estimation
- `RNG_SEED_MAGIC -666` — Numerical Recipes convention; triggers re-seeding
- `MAX_SIGMA_CLIP_RETRIES 5` — iteration cap to prevent divergence
- `KERNEL_BASE_ORDER_FORMULA` — document the polynomial basis formula with Alard & Lupton reference

**D. Documenting Undocumented Algorithms**

For functions with complex multi-step algorithms (e.g., `getStampStats3`), add Doxygen
`@brief` and `@details` blocks that:

- Explain the algorithm in prose (not just parameter lists)
- Reference the paper equation or section (e.g., "Alard & Lupton 1998, Eq. 3")
- Note any non-obvious invariants (e.g., "assumes stamps are background-subtracted")
- Describe the purpose of flags and special cases (e.g., when renormalization is triggered)

Priority: Functions larger than 50 lines or those implementing algorithms not obvious from code.

**E. Standardizing Naming Conventions**

Document a single preferred convention in `CONTRIBUTING.md`:

- **Legacy code (globals.h, existing camelCase functions):** Accept as-is; incremental improvement acceptable
- **New code:** Use `snake_case` for functions and local variables; reserve camelCase for type names
- **Future refactoring:** Gradual migration to snake_case is preferred, but not required for small fixes

### Contributor Checklist

When writing or modifying code in HOTPANTS, ensure:

- Global variables follow the `t*/i*/m*` naming convention or are documented in a prefix legend
- Loop variables in nesting depth > 2 use descriptive names (not i, j, k)
- All magic numbers are replaced with `#define` constants in defaults.h with explanatory comments
- All functions longer than 50 lines have a Doxygen `@brief` and `@details` block
- All undocumented algorithms include a reference to Alard & Lupton (1998) or the relevant paper
- Parity checks, bit-twiddling logic, and other non-obvious operations have explanatory comments
- Functions are decomposed if nesting depth exceeds 3 or cyclomatic complexity exceeds 5
- Dead code is removed; debug output uses consistent logging macros with explicit verbose-level checks
- Code follows established naming conventions (camelCase for legacy, snake_case for new)
- Optimization-critical sections are commented with their expected performance profile

### Why Readability Matters for HOTPANTS

Readability improvements directly support the modernization goals:

1. **Python API Clarity:** The C core must be self-documenting for Doxygen + Sphinx integration; vague variable names and hidden algorithms hinder API documentation quality.
2. **Faster Onboarding:** Explicit naming and documented algorithms reduce the time new contributors spend deciphering code intent.
3. **Debugging and Optimization:** Clear algorithm documentation makes it easier to identify bottlenecks and verify correctness during performance profiling.
4. **Long-term Maintainability:** Well-named code and comprehensive comments reduce technical debt and prevent subtle bugs from creeping in during refactoring.
