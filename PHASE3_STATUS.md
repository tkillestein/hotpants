# Phase 3 Integration Testing — Status & Next Steps

## Current State (May 6, 2026)

### ✅ Completed

**Infrastructure:**
- Phase 1: ctypes bindings complete (33/33 unit tests pass)
- Phase 2: Build system documentation (CONTRIBUTING.md, README, pyproject.toml)
- Pydantic models for configuration with full validation
- Test fixtures corrected (star_field, make_image, run_hotpants signatures)
- Mock C library properly set up for unit tests

**Python API:**
- `KernelConfig`, `RegionLayout`, `NoiseThresholds` — fully validated
- `fit_kernel()`, `spatial_convolve()` — signatures correct, orchestration logic in place
- `KernelSolution` — result dataclass with proper serialization
- Error handling and input validation working correctly

### ⧉ In Progress

**C Function Integration:**
- Skeleton implementations of `build_stamps()`, `fit_kernel_c()` exist
- These don't yet call actual C functions for region iteration and fitting
- Integration tests fixed but skipped until C functions are implemented

## What's Needed to Complete Phase 3

### 1. Proper buildStamps() Wrapper (2–3 hours)

**Current state:** Allocates stamp array but doesn't build stamps

**What's needed:**
- Iterate over regions (loop: nRegX × nRegY)
- For each region:
  - Define region bounds (rXMin, rXMax, rYMin, rYMax)
  - Call C `buildStamps()` with parameters matching main.c usage
  - Call C `getPsfCenters()` to locate bright stars
  - Handle memory for auxiliary stamp arrays
- Return complete stamp array ready for fitting

**Reference:** `src/main.c` lines 900–1000 (region loop and buildStamps calls)

**Signature to match:**
```c
void buildStamps(int sXMin, int sXMax, int sYMin, int sYMax, int* rPixX,
                 int* rPixY, int nStampX, int nStampY, int hwKSStamp,
                 stamp_struct* stamps, stamp_struct* tStamps,
                 stamp_struct* iStamps, float* iData, float* tData,
                 float tUThresh, float tLThresh);
```

### 2. Proper fitKernel_c() Wrapper (2–3 hours)

**Current state:** Skeleton that doesn't call C fitKernel()

**What's needed:**
- Call C `fitKernel()` with:
  - stamp array (from buildStamps)
  - template and science image pointers
  - optional noise image
  - output buffers for kernel coefficients
- Extract results:
  - kernel coefficients (size: nCompKer × nComp)
  - mean significance, scatter, chi2
  - number of skipped stamps
- Proper error handling

**Reference:** `src/main.c` lines 1150–1160 (fitKernel call and result handling)

**Signature to match:**
```c
void fitKernel(stamp_struct* stamps, float* imRef, float* imConv,
               float* imNoise, double* kernel_coeffs, double* meansig,
               double* scatter, int* n_skipped);
```

### 3. Enable & Run Integration Tests (1–2 hours)

**Once C functions are implemented:**
1. Remove `@pytest.mark.skip` from TestPythonAPIvsCLI
2. Run integration tests:
   ```bash
   export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH
   uv run pytest tests/test_api_integration.py::TestPythonAPIvsCLI -v
   ```
3. Debug failures and iterate
4. Compare Python API vs CLI outputs (tolerance: rtol=1e-4, atol=1e-3)

**Expected behavior:**
- Identical images → difference is ~noise (chi2 ≈ 0)
- Broadened PSF → kernel_norm > 1
- Narrowed PSF → kernel_norm < 1
- Python API and CLI produce numerically similar results

### 4. Add Regression Tests (1 hour)

**Current:**  `tests/test_regression.py` exists but is empty

**Implement:**
- Small synthetic image test case
- Run Python API → capture outputs
- Run C CLI → capture outputs  
- Assert outputs match within tolerance
- Document expected values as regression baseline

## Implementation Checklist

- [ ] Study `src/main.c` region loop (lines 900–1000)
- [ ] Study `src/alard.c` for buildStamps/fitKernel details
- [ ] Implement region iteration in `build_stamps()`
- [ ] Implement buildStamps C call with proper parameters
- [ ] Implement fitKernel C call with output handling
- [ ] Test with single region (nRegX=1, nRegY=1)
- [ ] Test with multiple regions
- [ ] Unskip integration tests
- [ ] Debug and fix test failures
- [ ] Validate numerical accuracy vs CLI
- [ ] Implement regression tests
- [ ] Document any limitations or assumptions

## Key Files

**To modify:**
- `src/hotpants/_core.py` — `build_stamps()` and `fit_kernel_c()` functions

**To read:**
- `src/main.c` — lines 900–1000 (region iteration, buildStamps pattern)
- `src/alard.c` — fitKernel implementation
- `src/functions.c` — buildStamps, getPsfCenters implementation
- `hotpants_api.h` — C function signatures

**Tests:**
- `tests/test_api_integration.py` — integration tests (fixed, ready to run)
- `tests/test_regression.py` — regression tests (empty, needs implementation)

## Estimated Effort

- **Analysis & Setup:** 1 hour (understand C flow)
- **Implementation:** 5 hours (buildStamps + fitKernel)
- **Testing & Debugging:** 2 hours (enable tests, fix failures)
- **Regression Tests:** 1 hour
- **Documentation:** 0.5 hours
- **Total: 9.5 hours**

## Critical Success Factors

1. **Memory management:** Proper allocation/freeing of stamp structures
2. **Region iteration:** Correct handling of multiple regions
3. **Parameter mapping:** Global variables → C function parameters
4. **Output handling:** Extraction of results from C via pointers
5. **Numerical accuracy:** Tolerance matching for floating-point operations

## Notes for Future Implementation

1. Start with single-region tests (nRegX=1, nRegY=1)
2. Use small images (128×128) for fast debugging
3. Add logging/debug output for parameter inspection
4. Compare intermediate results vs main.c at each step
5. Profile performance once functional

---

**Status:** Phase 3.1 (fixtures) complete. Phase 3.2 (C integration) blocked pending implementation. Phase 3.3–3.4 (tests) ready to run once C functions work.

The Python API is **production-ready from a structural/validation perspective**. Full end-to-end testing requires completing the C function orchestration.
