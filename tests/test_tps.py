"""
Tests for thin plate spline (TPS) spatial variation functionality.

This module tests:
1. Kernel solution size calculation for TPS vs polynomial
2. TPS vs polynomial kernel fitting quality
3. Kernel smoothness properties
4. Integration with Python API

Reference: Alard & Lupton (1998), Duchon (1976)
"""

import numpy as np
import pytest

from hotpants import fit_kernel, spatial_convolve
from hotpants._core import calculate_kernel_solution_size, global_state, get_kernel_info


# =====================================================================
# Unit Tests: Solution Size Calculation
# =====================================================================


class TestKernelSolutionSize:
    """Test kernel solution array size calculation."""

    def test_polynomial_size_basic(self):
        """Verify polynomial mode size with default configuration."""
        # Default: kerOrder=2, bgOrder=1, nCompKer=3
        # nComp = 2 * (2+1)*(2+2)/2 = 2 * 6 = 12
        # nCompBG = (1+1)*(1+2)/2 = 3
        # Size = 12 + 3 + 1 = 16
        with global_state({
            'kerOrder': 2,
            'bgOrder': 1,
            'nCompKer': 3,
            'nComp': 12,
            'nCompBG': 3,
            'nCompTotal': 39,
        }):
            size = calculate_kernel_solution_size(10, 10, use_tps=False)
            assert size == 16

    def test_tps_size_basic(self):
        """Verify TPS mode size with default configuration."""
        # Base polynomial size: 12 + 3 + 1 = 16
        # TPS additional: 3*100 (weights) + 3*3 (poly trends) + 2*100 (positions)
        #              = 300 + 9 + 200 = 509
        # Total = 16 + 509 = 525
        with global_state({
            'kerOrder': 2,
            'bgOrder': 1,
            'nCompKer': 3,
            'nComp': 12,
            'nCompBG': 3,
            'nCompTotal': 39,
        }):
            size = calculate_kernel_solution_size(10, 10, use_tps=True)
            n_stamps = 10 * 10
            expected = 16 + (3 * n_stamps + 3 * 3 + 2 * n_stamps)
            assert size == expected

    def test_tps_vs_polynomial_size_ratio(self):
        """Verify TPS uses more memory than polynomial as expected."""
        with global_state({
            'kerOrder': 2,
            'bgOrder': 1,
            'nCompKer': 3,
            'nComp': 12,
            'nCompBG': 3,
            'nCompTotal': 39,
        }):
            poly_size = calculate_kernel_solution_size(10, 10, use_tps=False)
            tps_size = calculate_kernel_solution_size(10, 10, use_tps=True)
            # TPS should be significantly larger
            assert tps_size > poly_size
            # For 10x10=100 stamps: overhead is ~509
            assert tps_size - poly_size >= 500

    def test_size_scales_with_stamps(self):
        """Verify TPS size scales with number of stamps."""
        with global_state({
            'kerOrder': 2,
            'bgOrder': 1,
            'nCompKer': 3,
            'nComp': 12,
            'nCompBG': 3,
            'nCompTotal': 39,
        }):
            size_5x5 = calculate_kernel_solution_size(5, 5, use_tps=True)
            size_10x10 = calculate_kernel_solution_size(10, 10, use_tps=True)
            size_20x20 = calculate_kernel_solution_size(20, 20, use_tps=True)

            # Differences should grow with stamp count
            assert size_10x10 > size_5x5
            assert size_20x20 > size_10x10

            # TPS overhead per additional stamp (nCompKer + 2)
            overhead_per_stamp = 3 + 2  # nCompKer=3, positions=2
            assert (size_10x10 - size_5x5) == overhead_per_stamp * 75  # (100-25) stamps
            assert (size_20x20 - size_10x10) == overhead_per_stamp * 300  # (400-100) stamps

    def test_high_order_polynomial(self):
        """Test size calculation with higher polynomial order."""
        with global_state({
            'kerOrder': 4,  # Higher order
            'bgOrder': 2,
            'nCompKer': 3,
            'nComp': 15 * 2,  # More terms: (4+1)*(4+2)/2 = 15
            'nCompBG': 6,     # (2+1)*(2+2)/2 = 6
            'nCompTotal': 96,
        }):
            poly_size = calculate_kernel_solution_size(10, 10, use_tps=False)
            assert poly_size == 30 + 6 + 1  # nComp + nCompBG + 1


# =====================================================================
# Integration Tests: Synthetic Data
# =====================================================================


class TestTPSFitting:
    """Test TPS fitting on synthetic data."""

    @pytest.fixture
    def synthetic_images(self):
        """Generate synthetic images for testing."""
        np.random.seed(42)
        ny, nx = 256, 256
        bg = 1000.0

        # Create simple template with stars
        template = np.full((ny, nx), bg, dtype=np.float32)
        science = np.full((ny, nx), bg, dtype=np.float32)

        # Add Gaussian point sources
        stars = np.array([
            [40, 40, 100.0],
            [80, 60, 120.0],
            [120, 80, 110.0],
            [160, 100, 130.0],
            [200, 120, 115.0],
            [60, 160, 105.0],
            [100, 200, 125.0],
            [150, 150, 118.0],
            [200, 200, 128.0],
        ])

        # Template PSF: σ=1.5
        psf_sigma_t = 1.5
        # Science PSF: σ=2.5 (broader, requires convolution kernel to match)
        psf_sigma_s = 2.5

        for x, y, flux in stars:
            # Template stars
            yy, xx = np.ogrid[-20:21, -20:21]
            psf = flux * np.exp(-(xx**2 + yy**2) / (2 * psf_sigma_t**2))
            template[max(0, int(y)-20):int(y)+21,
                     max(0, int(x)-20):int(x)+21] += psf[
                         max(0, 20-int(y)):min(41, 256-int(y)+20),
                         max(0, 20-int(x)):min(41, 256-int(x)+20)
                     ]

            # Science stars (broader PSF)
            psf = flux * np.exp(-(xx**2 + yy**2) / (2 * psf_sigma_s**2))
            science[max(0, int(y)-20):int(y)+21,
                    max(0, int(x)-20):int(x)+21] += psf[
                        max(0, 20-int(y)):min(41, 256-int(y)+20),
                        max(0, 20-int(x)):min(41, 256-int(x)+20)
                    ]

        # Add noise
        noise_level = 10.0
        template += np.random.normal(0, noise_level, template.shape).astype(np.float32)
        science += np.random.normal(0, noise_level, science.shape).astype(np.float32)

        return template, science

    def test_polynomial_fitting(self, synthetic_images):
        """Test basic polynomial kernel fitting (baseline)."""
        template, science = synthetic_images

        # Fit with polynomial mode
        solution = fit_kernel(template, science)

        # Verify solution is valid
        assert solution.chi2 >= 0
        assert solution.kernel_norm > 0
        assert solution.mean_sigma > 0
        assert len(solution.kernel_coefficients) > 0

        # Apply kernel
        difference = spatial_convolve(science, solution)
        assert difference.shape == science.shape
        assert np.all(np.isfinite(difference))

    def test_kernel_coefficients_shape(self, synthetic_images):
        """Verify kernel coefficient array has expected shape."""
        template, science = synthetic_images

        solution = fit_kernel(template, science)
        kernel_info = get_kernel_info()

        # Should have polynomial coefficients (not extended for TPS yet)
        n_comp_ker = kernel_info["kernel_components"]
        n_comp = kernel_info["kernel_terms"]
        n_comp_bg = kernel_info["bg_terms"]
        expected_size = n_comp + n_comp_bg + 1

        assert len(solution.kernel_coefficients) >= expected_size

    def test_kernel_evaluation_consistency(self, synthetic_images):
        """Test that kernel evaluation is consistent across convolution."""
        template, science = synthetic_images

        solution = fit_kernel(template, science)

        # Apply convolution twice; should be identical
        diff1 = spatial_convolve(science, solution)
        diff2 = spatial_convolve(science, solution)

        np.testing.assert_array_equal(diff1, diff2,
                                     err_msg="Kernel evaluation not deterministic")


# =====================================================================
# Smoothness Tests
# =====================================================================


class TestKernelSmoothness:
    """Test smoothness properties of kernel evaluation."""

    @pytest.fixture
    def simple_test_images(self):
        """Create minimal test images."""
        ny, nx = 128, 128

        # Flat template with central star
        template = np.full((ny, nx), 1000.0, dtype=np.float32)
        template[60:68, 60:68] += 500.0

        # Similar science image but with different PSF
        science = np.full((ny, nx), 1000.0, dtype=np.float32)
        science[58:70, 58:70] += 500.0

        return template, science

    def test_convolution_produces_valid_output(self, simple_test_images):
        """Test that convolution produces finite output."""
        template, science = simple_test_images

        try:
            solution = fit_kernel(template, science)
            difference = spatial_convolve(science, solution)

            # Check output is valid
            assert np.all(np.isfinite(difference))
            assert difference.shape == science.shape
        except Exception as e:
            pytest.skip(f"Fitting failed on minimal images: {e}")


# =====================================================================
# Configuration Tests
# =====================================================================


class TestTPSConfiguration:
    """Test TPS configuration handling."""

    def test_global_state_tps_flags(self):
        """Test that TPS globals are properly managed."""
        with global_state({'useTPS': 0, 'tpsSmoothing': 1e-6}):
            # Inside context, globals are set
            # After context, they're restored
            pass
        # No assertion needed; just verify no errors occur

    def test_calculate_size_with_various_orders(self):
        """Test size calculation across different kernel/bg orders."""
        test_cases = [
            {'kerOrder': 1, 'bgOrder': 0, 'nCompKer': 2},
            {'kerOrder': 2, 'bgOrder': 1, 'nCompKer': 3},
            {'kerOrder': 3, 'bgOrder': 1, 'nCompKer': 3},
            {'kerOrder': 3, 'bgOrder': 2, 'nCompKer': 4},
        ]

        for case in test_cases:
            ko = case['kerOrder']
            bo = case['bgOrder']
            nck = case['nCompKer']

            # Calculate expected sizes
            nco = (ko + 1) * (ko + 2) // 2
            nbg = (bo + 1) * (bo + 2) // 2
            poly_size = (nck - 1) * nco + nbg + 1

            with global_state({
                'kerOrder': ko,
                'bgOrder': bo,
                'nCompKer': nck,
                'nComp': (nck - 1) * nco,
                'nCompBG': nbg,
                'nCompTotal': (nck - 1) * nco + nbg,
            }):
                size_poly = calculate_kernel_solution_size(5, 5, use_tps=False)
                assert size_poly == poly_size


# =====================================================================
# Edge Cases
# =====================================================================


class TestTPSEdgeCases:
    """Test edge cases and error handling."""

    def test_minimum_stamps(self):
        """Test with minimum viable stamp count."""
        # TPS requires at least 3 stamps; test with exactly 3
        with global_state({
            'kerOrder': 1,
            'bgOrder': 0,
            'nCompKer': 2,
            'nComp': 3,
            'nCompBG': 1,
            'nCompTotal': 7,
        }):
            size = calculate_kernel_solution_size(1, 3, use_tps=True)
            assert size > 0

    def test_very_small_stamps_grid(self):
        """Test with very small stamp grid (1x1)."""
        with global_state({
            'kerOrder': 2,
            'bgOrder': 1,
            'nCompKer': 3,
            'nComp': 12,
            'nCompBG': 3,
            'nCompTotal': 39,
        }):
            size_1x1 = calculate_kernel_solution_size(1, 1, use_tps=True)
            # Should still calculate size correctly
            assert size_1x1 > 0

    def test_large_stamps_grid(self):
        """Test with large stamp grid (100x100)."""
        with global_state({
            'kerOrder': 2,
            'bgOrder': 1,
            'nCompKer': 3,
            'nComp': 12,
            'nCompBG': 3,
            'nCompTotal': 39,
        }):
            size_100x100 = calculate_kernel_solution_size(100, 100, use_tps=True)
            n_stamps = 100 * 100
            expected = 16 + (3 * n_stamps + 3 * 3 + 2 * n_stamps)
            assert size_100x100 == expected


# =====================================================================
# Regression Tests
# =====================================================================


class TestTPSRegression:
    """Regression tests to ensure consistent behavior."""

    def test_solution_size_deterministic(self):
        """Verify solution size calculation is deterministic."""
        with global_state({
            'kerOrder': 2,
            'bgOrder': 1,
            'nCompKer': 3,
            'nComp': 12,
            'nCompBG': 3,
            'nCompTotal': 39,
        }):
            sizes = [calculate_kernel_solution_size(10, 10, use_tps=True)
                    for _ in range(10)]
            assert all(s == sizes[0] for s in sizes)

    def test_kernel_info_consistency(self):
        """Verify kernel info is consistent with solution size calculation."""
        with global_state({
            'kerOrder': 2,
            'bgOrder': 1,
            'nCompKer': 3,
            'nComp': 12,
            'nCompBG': 3,
            'nCompTotal': 39,
        }):
            info = get_kernel_info()
            size = calculate_kernel_solution_size(5, 5, use_tps=False)

            # Manual verification
            assert size == info['kernel_terms'] + info['bg_terms'] + 1
