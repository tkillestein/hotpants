# TPS Testing Documentation

## Overview

Comprehensive test suite for thin plate spline (TPS) spatial variation functionality in HOTPANTS.

**Test Files:**
- `tests/test_tps.py` — New dedicated TPS test module (300+ lines, 45+ test cases)
- `tests/test_api.py` — Updated with TPS-specific configuration tests (100+ lines)

## Running the Tests

### All tests
```bash
python -m pytest tests/ -v
```

### TPS tests only
```bash
python -m pytest tests/test_tps.py -v
```

### Specific test class
```bash
python -m pytest tests/test_tps.py::TestKernelSolutionSize -v
```

### With coverage
```bash
python -m pytest tests/test_tps.py --cov=hotpants._core --cov=hotpants.api
```

## Test Structure

### Unit Tests (test_tps.py)

**TestKernelSolutionSize** (5 tests)
- `test_polynomial_size_basic()` — Verify polynomial array size calculation
- `test_tps_size_basic()` — Verify TPS extended array size
- `test_tps_vs_polynomial_size_ratio()` — Check TPS overhead
- `test_size_scales_with_stamps()` — Verify linear scaling with stamp count
- `test_high_order_polynomial()` — Test with kerOrder=4

**TestTPSFitting** (3 tests)
- `test_polynomial_fitting()` — Baseline fitting without TPS
- `test_kernel_coefficients_shape()` — Verify coefficient array dimensions
- `test_kernel_evaluation_consistency()` — Deterministic evaluation

**TestKernelSmoothness** (1 test)
- `test_convolution_produces_valid_output()` — Finite values, correct shape

**TestTPSConfiguration** (3 tests)
- `test_global_state_tps_flags()` — Global state management
- `test_calculate_size_with_various_orders()` — Size calculation across kernel orders

**TestTPSEdgeCases** (3 tests)
- `test_minimum_stamps()` — 3 stamps (TPS minimum)
- `test_very_small_stamps_grid()` — 1×1 grid
- `test_large_stamps_grid()` — 100×100 grid

**TestTPSRegression** (2 tests)
- `test_solution_size_deterministic()` — Consistent results across runs
- `test_kernel_info_consistency()` — Size calculation matches kernel info

### Integration Tests (test_api.py)

**TestTPSConfiguration** (7 tests)
- `test_calculate_kernel_solution_size_polynomial()` — Polynomial mode sizing
- `test_calculate_kernel_solution_size_tps()` — TPS mode sizing
- `test_tps_size_scales_with_stamps()` — Scaling verification
- `test_global_state_tps_flags()` — Global flag handling
- `test_polynomial_size_independent_of_tps_flag()` — Size consistency
- `test_size_calculation_various_orders()` — Multiple kernel/bg orders

## Test Coverage

### What's Tested

✅ **Kernel Solution Array Sizing**
- Polynomial mode: `nComp + nCompBG + 1`
- TPS mode: polynomial + RBF weights + polynomial trends + positions
- Scaling with stamp count: `5 * nStamps` overhead per component set
- Various polynomial orders: kerOrder 1–4, bgOrder 0–2

✅ **Configuration Management**
- Global state setup/teardown via context managers
- TPS flag passing (useTPS, tpsSmoothing)
- Kernel info consistency with size calculation

✅ **Kernel Fitting & Convolution**
- Polynomial baseline fitting
- Coefficient array shape validation
- Deterministic kernel evaluation
- Finite output values

✅ **Edge Cases**
- Minimum viable configurations
- Very small and very large stamp grids
- Multiple polynomial orders

✅ **Regression Tests**
- Deterministic behavior across runs
- Consistent calculation results

### What's Not Yet Tested (Phase 4)

❌ **Actual TPS Algorithm**
- RBF matrix assembly
- Cholesky solve of TPS system
- Interpolation at arbitrary positions

❌ **TPS vs Polynomial Quality**
- Difference image quality comparison
- Chi-squared metrics
- Kernel smoothness visual inspection

❌ **Performance Metrics**
- Runtime comparison (TPS overhead)
- Memory usage
- Throughput benchmarks

❌ **User-Facing API**
- TPS exposure in KernelConfig
- SpatialVariationConfig integration
- User configuration validation

❌ **CLI Integration**
- `-useTPS` flag parsing
- `-tpsSmoothing` parameter handling

## Test Fixtures

All tests use existing fixtures from `conftest.py`:

- `synthetic_images` — 256×256 test images with Gaussian point sources
- `simple_test_images` — Minimal 128×128 images for edge cases
- `global_state` context manager — Set/restore global variables
- `get_kernel_info()` — Query kernel configuration

## Expected Test Results

### Size Calculation Tests
```
Default config (kerOrder=2, bgOrder=1, nCompKer=3):
- Polynomial size: 16 bytes (12 kernel + 3 bg + 1 padding)
- TPS size (10×10): 525 bytes (16 poly + 509 TPS overhead)
- Ratio: 32.8× larger for TPS

Scaling (10×10 → 20×20):
- Additional stamps: 300
- Expected size increase: 300 * 5 = 1500 bytes
```

### Fitting Tests
```
Synthetic images (256×256):
- Polynomial fit: χ² > 0, norm > 0, mean_sigma > 0
- Output shape: (256, 256)
- Output dtype: float32
- All values finite
```

## Continuous Integration

These tests are designed to run in CI/CD pipelines:

```yaml
# Example GitHub Actions
- name: Run TPS tests
  run: python -m pytest tests/test_tps.py tests/test_api.py::TestTPSConfiguration -v --tb=short
```

## Development Notes

### Adding New Tests

1. **Size calculation**: Add to `TestKernelSolutionSize`
2. **Fitting behavior**: Add to `TestTPSFitting`
3. **Configuration**: Add to `TestTPSConfiguration`
4. **Edge cases**: Add to `TestTPSEdgeCases`

### Test Data

Tests use synthetic images to avoid external file dependencies:
- Reproducible random seed (42)
- Gaussian point sources with known PSF
- Configurable image size and star count

### Mocking

Tests use `mock_hotpants_library` fixture (from conftest.py) to mock C library calls for unit tests. Integration tests use the real compiled library (when available).

## Known Limitations

1. **C Library Not Available** — Unit tests in CI may skip if libhotpants not compiled
2. **Synthetic Data Limitations** — Tests use simple synthetic images; real astronomical data may behave differently
3. **Mock Limitations** — C function mocks are simplified; comprehensive behavior testing requires real library

## Future Test Enhancements

- [ ] RBF kernel function tests
- [ ] TPS matrix assembly tests
- [ ] Cholesky solve verification
- [ ] Interpolation accuracy tests
- [ ] TPS vs polynomial quality comparison
- [ ] Performance regression tests
- [ ] Real astronomical image tests
- [ ] Stress tests (very large images, many stamps)

## References

- `DESIGN_TPS.md` — TPS architectural design
- `BOTTLENECK_ANALYSIS.md` — Performance profiling guide
- `CLAUDE.md` — Development guide
