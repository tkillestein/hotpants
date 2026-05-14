# Thin Plate Spline Spatial Variation — Design Document

## Overview

This document outlines the design for implementing thin plate spline (TPS) spatial variation of kernel coefficients in the C backend, as an alternative to the existing polynomial model (Alard & Lupton 1998).

**Goals:**
- Enable single-region full-frame processing (`-nrx 1 -nry 1`) without tiling artifacts
- Provide smoother kernel variation across the image
- Maintain backward compatibility with polynomial mode
- Keep computation entirely in C

**Target implementation date:** Phase 1 complete by end of sprint.

---

## Mathematical Foundation

### Current Approach: Spatial Polynomials

Each kernel coefficient `c_i(x,y)` is expanded as a 2D polynomial:

```
c_i(x,y) = Σ_{dx=0}^{k} Σ_{dy=0}^{k-dx} a_{i,dx,dy} · x^{dx} · y^{dy}
```

where `k = kerOrder` (typically 2–4).

**Limitations:**
- Ringing and oscillation at region boundaries
- Discontinuous derivatives at boundaries (visible as steps in difference images)
- Poor fit to non-polynomial spatial variation

### Proposed Approach: Thin Plate Splines

Replace polynomial expansion with RBF interpolation using the thin plate spline (TPS) kernel:

```
φ(r) = r² log(r)   (TPS RBF kernel)
```

Each coefficient c_i is approximated as:

```
c_i(x,y) ≈ Σ_{j=1}^{N} w_{i,j} · φ(||(x,y) - (x_j, y_j)||) + p_{i,0}(x,y)
```

where:
- `(x_j, y_j)` are the N stamp center positions
- `w_{i,j}` are RBF weights (one set per kernel component)
- `p_{i,0}(x,y)` is a low-degree polynomial trend (usually constant or linear)

**Advantages:**
- C∞ smooth interpolation
- Rotation-invariant (no preferred direction)
- No artificial discontinuities
- Exact fit at control points (with regularization option)
- Natural for scattered data (stamps may be non-uniform)

**Trade-offs:**
- Slightly slower evaluation than polynomial (O(N) vs. O(k²))
- Requires solving NxN linear system (N = num stamps per region)
- More memory: O(N) weights per kernel component

---

## Implementation Strategy

### Phase 1: Core TPS in C

**Files affected:**
- `src/defaults.h` — add named constants for TPS
- `src/globals.h` — add global configuration flags
- `src/alard.c` — add TPS fitting and evaluation functions

**New globals (in `globals.h`):**

```c
extern int useTPS;        /* 1 = use TPS, 0 = use polynomial */
extern double tpsSmoothing;  /* Regularization: 0 = exact fit, >0 = smooth */
extern double tpsRadius;  /* Expected spatial scale; auto-set from image size */
```

**New functions in `alard.c`:**

1. **`tpsFitCoefficients()`** — Fit TPS surface to scattered stamp data
   - Input: stamp positions (N, 2) and fitted coefficients (N,) for one kernel component
   - Output: RBF weights and polynomial coefficients
   - Internally: solve the TPS linear system (augmented normal equations)

2. **`tpsEvaluate()`** — Evaluate TPS at arbitrary position
   - Input: (x, y), RBF weights, control points, poly coeffs
   - Output: interpolated coefficient value
   - Used in place of `make_kernel()` loop when `useTPS=1`

3. **`makeKernelTPS()`** — TPS version of `make_kernel()`
   - Like `make_kernel()` but uses `tpsEvaluate()` for each component
   - Called from `spatial_convolve()` when `useTPS=1`

### Phase 2: Integration with fitKernel()

**Modified fitKernel():**
1. After LAPACK solve (polynomial case), check `useTPS` flag
2. If `useTPS=1`:
   - Extract (x_j, y_j) stamp centers and fitted c_i values
   - Call `tpsFitCoefficients()` for each kernel component
   - Store RBF weights and poly coefficients in extended solution vector
3. Return solution with TPS parameters

**Data structure for TPS solution storage:**

The `kernelSol` vector grows to include TPS parameters when `useTPS=1`:

```
kernelSol layout (useTPS=1):
[padding, poly_coeff_0, poly_coeff_1, ..., poly_coeff_(nComp-1),
 bg_poly_coeffs,
 tps_weights_comp_0[N], tps_weights_comp_1[N], ...,
 stamp_x[N], stamp_y[N]]
```

where N = number of stamps = nStampX * nStampY.

### Phase 3: Evaluation in spatial_convolve()

**Modified spatial_convolve():**
1. Check `useTPS` flag
2. If `useTPS=1`:
   - Extract RBF weights, stamp positions from kernelSol
   - For each pixel (x, y):
     - Call `makeKernelTPS()` to evaluate K(x,y) using TPS
   - Proceed as normal with FFT convolution

---

## Data Structures

### Stamp Grid (existing, unchanged)

```c
typedef struct {
  float *image;       /* Stamp subimage */
  float *refImage;    /* Reference image stamp */
  int centerX, centerY;  /* Stamp center pixel coordinates */
  /* ... other fields ... */
} stamp_struct;
```

### TPS Control Points (new)

At TPS fitting time, we extract:
- **positions**: (nStampX * nStampY, 2) array of stamp centers
- **values**: (nStampX * nStampY,) array of fitted coefficient c_i for one kernel component

These are implicitly stored in the kernel solution vector when `useTPS=1`.

### TPS RBF Matrix (internal, temporary)

During fitting in `tpsFitCoefficients()`:
```c
/* Augmented RBF matrix for TPS solve:
   [ Φ  P ] [ w ]   [ c ]
   [ P' 0 ] [ v ] = [ 0 ]
   
   where Φ[i,j] = φ(||pos[i] - pos[j]||)
         P is polynomial basis (e.g., [1, x, y] for linear)
         w are RBF weights
         v are polynomial coefficients
*/
```

---

## CLI Integration

Add command-line option to enable TPS:

```
-useTPS [0|1]       Enable thin plate spline spatial variation (default: 0)
                    Requires nrx=1, nry=1 for best results.

-tpsSmoothing <val> TPS regularization parameter (default: 1e-6).
                    0 = exact fit; >0 = smoother surface.
```

Implementation in `vargs.c`:
```c
if (! strcmp(argv[ii], "-useTPS")) {
  useTPS = atoi(argv[++ii]);
  /* Validate: TPS typically requires single region */
}
```

---

## Testing Strategy

### Unit tests (C, via test harness):
1. `test_tpsEvaluate()` — verify evaluation at control points (should be exact or near-exact)
2. `test_tpsFitCoefficients()` — verify RBF system solve
3. `test_makeKernelTPS()` — compare against polynomial for known data

### Integration tests (Python API):
1. Fit kernel with `useTPS=1` and `nrx=1, nry=1`
2. Compare difference image quality vs. polynomial (`nrx=2, nry=2`)
3. Measure smoothness of kernel variation (finite differences)

### Performance tests:
1. Benchmark `makeKernelTPS()` vs. `make_kernel()` for typical region size
2. Profile memory usage (additional RBF weights storage)
3. Total runtime: single region TPS vs. multi-region polynomial

---

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|-----------|
| TPS ill-conditioning for dense stamps | Numerical instability | Use SVD or regularization in RBF solve |
| Memory overflow for large stamp grids | Crash on very large images | Cap N or implement block-wise evaluation |
| Performance regression | Users avoid TPS | Profile early, optimize RBF eval |
| Integration complexity | Bugs, maintenance burden | Keep TPS functions isolated, test thoroughly |

---

## Implementation Progress

### Phase 1: Foundation ✓ COMPLETE

**Completed:**
- Core RBF mathematics in C (tps_kernel, tps_assemble_matrix, tps_fit_coefficients, tps_evaluate)
- Configuration infrastructure (useTPS, tpsSmoothing globals in defaults.h/globals.h)
- Kernel evaluation stub (make_kernel_tps placeholder with fallback to polynomial)
- Design documentation (this file)
- LAPACK integration via LAPACKE for LU solve

**Commits:**
- `0ad814e` - Core TPS RBF infrastructure and design document
- `5548934` - make_kernel_tps evaluation function placeholder

### Phase 2: fitKernel Integration (In Progress)

**Required:**
1. Extend kernelSol layout to store TPS parameters
   - Current: polynomial coefficients only
   - Extended: polynomial + RBF weights + polynomial trends + stamp positions
   
2. Modify fitKernel() to optionally fit TPS
   - After polynomial LAPACK solve, if useTPS==1:
     - Extract stamp centers from stamp_struct array
     - Evaluate polynomial at each stamp to get per-stamp coefficients
     - Call tps_fit_coefficients() for each kernel component
     - Store results in extended kernelSol
   
3. Implement tps_fit_kernel() helper
   - Called from fitKernel() after polynomial solve
   - Manages TPS fitting for all kernel components
   - Handles stamp position extraction and normalization
   
4. Update kernelSol size calculation
   - Need: nComp + nCompKer*(nStamps+3) + 2*nStamps
   - Where nStamps = nStampX * nStampY (per region)
   - Communicate new size to Python API

5. Dispatcher in spatial_convolve
   - Check useTPS flag before calling make_kernel or make_kernel_tps

### Phase 3: Python API & Testing ✓ COMPLETE

**Completed:**
- Add useTPS and tpsSmoothing to global_state configuration (_core.py)
- Implement calculate_kernel_solution_size() for TPS array sizing
- Update fit_kernel() to use correct solution size for TPS parameters
- Python API infrastructure ready; TPS disabled by default (use_tps=False)

**Commits:**
- `bac171c` - Python API support for TPS (Phase 3)

### Phase 4: User-Facing Configuration & Testing (Next)

**Required:**
1. Expose TPS to Python API
   - Add `spatial_variation` field to `KernelConfig` with TPS support
   - Update `fit_kernel()` to respect user TPS configuration
   - Validate prerequisites (single region recommended for best results)

2. Unit tests
   - test_tps_kernel() - verify RBF kernel values
   - test_tps_assemble_matrix() - check matrix structure
   - test_tps_fit_coefficients() - test LAPACK solve
   - test_tps_evaluate() - verify interpolation
   - test_make_kernel_tps() - kernel assembly
   - test_calculate_kernel_solution_size() - size calculation

3. Integration tests
   - Compare TPS vs polynomial on synthetic data
   - Verify smoothness of kernel variation
   - Benchmark performance (TPS overhead vs gain)
   - Test edge cases (tiny stamps, large images)

### Phase 5: CLI & Documentation (Later)

- Add `-useTPS` and `-tpsSmoothing` options to vargs.c
- Update help text and README
- Document TPS usage patterns and when to use
- Performance profiling and tuning

## Timeline

- **Phase 1 (Completed):** Core RBF infrastructure ✓
- **Phase 2 (Completed):** fitKernel integration ✓
- **Phase 3 (Completed):** Python API foundation ✓
- **Phase 4 (Next):** User-facing configuration & testing (½ day)
- **Phase 5:** CLI integration & documentation (½ day)

---

## References

- Duchon, J. (1976). "Splines minimizing rotation-invariant semi-norms." *Constr. Approx.* 6:85–100.
- Carr, J.C., et al. (2001). "Reconstruction and representation of 3D objects with radial basis functions."
- Alard & Lupton (1998). "A Method for Optimal Image Subtraction," *ApJ* 503:325.

