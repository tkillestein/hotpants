# Phase 3 Integration Testing — Status & Next Steps

## Current State (May 6, 2026 - Late Evening)

### ✅ Completed

**Infrastructure:**
- Phase 1: ctypes bindings complete (33/33 unit tests pass)
- Phase 2: Build system documentation (CONTRIBUTING.md, README, pyproject.toml)
- Pydantic models for configuration with full validation
- Test fixtures corrected (star_field, make_image, run_hotpants signatures)
- Mock C library properly set up for unit tests
- Fixed integration test mocking to only apply to unit tests (not integration/regression)
- Integration test CLI option fixes (-ko, -bgo instead of -k, -bg)
- Integration test output comparison logic updated

**Python API:**
- `KernelConfig`, `RegionLayout`, `NoiseThresholds` — fully validated
- `fit_kernel()`, `spatial_convolve()` — signatures correct, orchestration logic in place
- `KernelSolution` — result dataclass with proper serialization
- Error handling and input validation working correctly

### ⧉ In Progress

**C Function Integration (build_stamps wrapper):**
- `build_stamps()` wrapper designed with proper region iteration
- Calls actual C `buildStamps()` from functions.c with correct signature
- **Key discovery:** C buildStamps expects region-extracted image data, not full images
- Integration tests blocked: segfault when calling buildStamps with full image pointers
- **Next step:** Implement region extraction before calling C function

## Critical Discovery: C buildStamps() Requires Region Data Extraction

The C `buildStamps()` function in functions.c (line 250) has the signature:
```c
void buildStamps(int sXMin, int sXMax, int sYMin, int sYMax, int* niS, int* ntS,
                 int getCenters, int rXBMin, int rYBMin, stamp_struct* ciStamps,
                 stamp_struct* ctStamps, float* iRData, float* tRData,
                 float hardX, float hardY);
```

**Critical parameter:** `rPixX` and `rPixY` are NOT pointers (as described in hotpants_api.h), but integers representing the region width and height. The function internally uses these to compute row-major offsets into `iRData` and `tRData`.

**This means:**
- `iRData` and `tRData` must point to **region-extracted image buffers**, not the full image
- The wrapper must call `cutRegion()` to extract regions before calling `buildStamps()`
- Current implementation passes full image pointers → causes segfault

**Example from main.c (lines 746-748):**
```c
tRData = (float*)calloc(rPixX * rPixY, sizeof(float));
iRData = (float*)calloc(rPixX * rPixY, sizeof(float));
...
cutRegion(tData, tRData, ...); // Extract region from full image
cutRegion(iData, iRData, ...);
```

---

## What's Needed to Complete Phase 3

### 1. Proper buildStamps() Wrapper (3–4 hours)

**Current state:** Wrapper structure in place, but segfaults because it passes full image data instead of region-extracted data.

**What's needed:**
1. **Region extraction:** Before calling C `buildStamps()`, extract region data:
   ```python
   # For each region:
   region_t_data = allocate_array((r_pix_y, r_pix_x), dtype=np.float32)
   region_i_data = allocate_array((r_pix_y, r_pix_x), dtype=np.float32)
   cutRegion(template, region_t_data, ...)  # Extract via C or numpy
   cutRegion(science, region_i_data, ...)
   ```

2. **Iterate over regions and stamps:**
   - For each region: extract region data from full image
   - For each stamp within region:
     - Calculate stamp bounds (relative to region boundaries)
     - Call C `buildStamps()` with **region data** (not full image)
     - Track stamp counters (niS, ntS) — only increment if substamps found

3. **Call actual C buildStamps signature:**
   ```c
   buildStamps(sXMin, sXMax, sYMin, sYMax, &niS, &ntS, getCenters,
               rXBMin, rYBMin, ciStamps, ctStamps, iRData, tRData,
               hardX, hardY);
   ```

**Key implementation details:**
- `rPixX`, `rPixY`: NOT pointers; must pass region dimensions as integers via global state
- `iRData`, `tRData`: Must point to region-extracted buffers
- Stamp counters (`niS`, `ntS`): Passed as pointers, incremented by C function
- `getCenters`: Pass 1 to auto-find bright centers (matching main.c default)

**Reference:** 
- `src/main.c` lines 746–748 (region buffer allocation)
- `src/main.c` lines 919–985 (region loop and per-stamp buildStamps calls)
- `src/functions.c` line 250 (actual buildStamps signature)
- `src/functions.c` line 443 (cutStamp function showing rPixX usage)
- `src/functions.c` line 99 (cutRegion function for extracting regions)

**Helper functions needed:**
- `cutRegion()` wrapper: Extract rectangular region from full image
  - Takes: full image, region bounds, output buffer
  - Returns: void (fills buffer in-place)
  - Reference: functions.c line 99

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

**Phase 3.2a — Region Data Extraction (NEW TASK, 2–3 hours):**
- [ ] Add `cutRegion()` ctypes wrapper to `_core.py`
- [ ] Update `build_stamps()` to allocate per-region buffers
- [ ] Extract regions from full images before calling C buildStamps
- [ ] Handle stamp counter arrays (niS, ntS) properly
- [ ] Test single-region case (nRegX=1, nRegY=1) on 256×256 image

**Phase 3.2b — Integration Test Fixes (1–2 hours):**
- [ ] Fix segfault in build_stamps() region data handling
- [ ] Unskip integration tests
- [ ] Run test_identical_images_api_vs_cli
- [ ] Debug numerical differences (tolerance: rtol=1e-4, atol=1e-3)
- [ ] Test broadened/narrowed PSF cases

**Phase 3.3 — fitKernel & Convolution (3–4 hours, if needed):**
- [ ] Study `src/alard.c` for fitKernel algorithm
- [ ] Implement `fit_kernel_c()` with proper output handling
- [ ] Implement `spatial_convolve_c()` (may be simpler, check if already working)
- [ ] Test with multi-region case (nRegX>1, nRegY>1)
- [ ] Validate numerical accuracy vs CLI

**Phase 3.4 — Regression & Documentation (1–2 hours):**
- [ ] Implement regression tests comparing API vs CLI output
- [ ] Document any limitations or assumptions
- [ ] Update PHASE3_STATUS.md with final status

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

## Effort Summary

**Completed (May 6, 2026):**
- ✓ Analysis & Setup: 2 hours (understand C flow, discover region data extraction requirement)
- ✓ Integration test fixtures fixed: 1 hour (CLI options, output comparison)
- ✓ build_stamps() structure: 2 hours (wrapper design, region iteration logic)

**Remaining (estimated):**
- Region data extraction: 2–3 hours (cutRegion wrapper, per-region buffers)
- buildStamps debugging: 2–3 hours (test single region, fix numerical issues)
- fit_kernel_c() & spatial_convolve_c(): 3–4 hours (if needed; spatial_convolve may work)
- Regression tests: 1–2 hours
- Documentation: 0.5 hours

**Total estimate: 11–15 hours** (original was 9.5; added complexity from region extraction discovery)

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
