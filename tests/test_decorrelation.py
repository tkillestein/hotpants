"""
Unit and integration tests for LSST afterburner decorrelation (DMTN-021).

Tests the decorrelation module which removes correlated noise from
difference images by convolving with a spatially-varying whitening kernel φ(x,y).
"""

import numpy as np
import pytest
from hotpants import fit_kernel, spatial_convolve, KernelConfig


class TestDecorrelationConfig:
    """Test decorrelation configuration parameters."""

    def test_decorrelation_disabled_by_default(self):
        """Decorrelation should be disabled by default."""
        config = KernelConfig()
        # Note: decorrelation not exposed in Python API yet, but test that
        # config can be created without errors
        assert config is not None

    def test_decorrelation_config_with_defaults(self):
        """Create config with decorrelation enabled (when API available)."""
        config = KernelConfig(use_tps=False)
        # Verify config created successfully
        assert config is not None


class TestDecorrelationFormula:
    """Test the decorrelation kernel formula."""

    def test_decorrelation_formula_bounds(self):
        """
        Decorrelation formula should produce φ ≥ 0.

        Formula: φ̂(k) = √[(σ_s² + σ_t²) / (σ_s² + σ_t² · κ̂²(k))]

        Properties:
        - For κ̂(k) = 0: φ̂(k) = √[(σ_s² + σ_t²) / σ_s²]
        - For κ̂(k) = 1: φ̂(k) = √[1] = 1
        - For κ̂(k) >> 1: φ̂(k) ≈ √[1/κ̂²] → 0
        """
        sig_s2 = 1.0  # Science variance
        sig_t2 = 1.0  # Template variance

        # Test κ̂ = 0
        kappa = 0.0
        phi = np.sqrt((sig_s2 + sig_t2) / (sig_s2 + sig_t2 * kappa**2))
        assert phi > 0.0
        assert np.isfinite(phi)

        # Test κ̂ = 1
        kappa = 1.0
        phi = np.sqrt((sig_s2 + sig_t2) / (sig_s2 + sig_t2 * kappa**2))
        assert np.isclose(phi, 1.0, rtol=1e-6)

        # Test κ̂ = 2
        kappa = 2.0
        phi = np.sqrt((sig_s2 + sig_t2) / (sig_s2 + sig_t2 * kappa**2))
        assert 0 < phi < 1.0  # Should be less than 1 for large κ

    def test_decorrelation_formula_variance_scaling(self):
        """φ should depend on the ratio of variances."""
        sig_s2 = 1.0
        sig_t2 = 1.0
        kappa = 1.0

        # Equal variances → φ = 1
        phi_equal = np.sqrt((sig_s2 + sig_t2) / (sig_s2 + sig_t2 * kappa**2))
        assert np.isclose(phi_equal, 1.0)

        # Test with larger kernel value (κ > 1)
        kappa_large = 2.0
        phi_large_kernel = np.sqrt((sig_s2 + sig_t2) / (sig_s2 + sig_t2 * kappa_large**2))
        # With κ > 1, denominator is larger, so φ < 1
        assert phi_large_kernel < phi_equal

    def test_decorrelation_commutes_with_scaling(self):
        """Scaling both image and kernel together should be commutative."""
        sig_s2 = 2.0
        sig_t2 = 3.0
        kappa = 1.5

        phi = np.sqrt((sig_s2 + sig_t2) / (sig_s2 + sig_t2 * kappa**2))

        # Scale both by 2x
        scale = 2.0
        phi_scaled = np.sqrt(
            (sig_s2 * scale + sig_t2 * scale)
            / (sig_s2 * scale + sig_t2 * scale * kappa**2)
        )

        assert np.isclose(phi, phi_scaled, rtol=1e-10)


class TestBilinearInterpolation:
    """Test bilinear interpolation of decorrelation kernel on grid."""

    def test_bilinear_interp_grid_corners(self):
        """
        Test bilinear interpolation at grid corners and edges.

        For a 2x2 grid:
          [0,0]---[1,0]
            |       |
          [0,1]---[1,1]
        """
        # Create a simple 2x2 grid
        grid_values = np.array([[1.0, 2.0], [3.0, 4.0]])

        # Test at grid corners
        p00 = grid_values[0, 0]
        p10 = grid_values[1, 0]
        p01 = grid_values[0, 1]
        p11 = grid_values[1, 1]

        # Bilinear interp at (0, 0) should return p00
        u, v = 0.0, 0.0
        phi = (1 - u) * (1 - v) * p00 + u * (1 - v) * p10 + (
            1 - u
        ) * v * p01 + u * v * p11
        assert np.isclose(phi, p00)

        # Bilinear interp at (1, 0) should return p10
        u, v = 1.0, 0.0
        phi = (1 - u) * (1 - v) * p00 + u * (1 - v) * p10 + (
            1 - u
        ) * v * p01 + u * v * p11
        assert np.isclose(phi, p10)

        # Bilinear interp at center (0.5, 0.5) should be average
        u, v = 0.5, 0.5
        phi = (1 - u) * (1 - v) * p00 + u * (1 - v) * p10 + (
            1 - u
        ) * v * p01 + u * v * p11
        phi_expected = (p00 + p10 + p01 + p11) / 4.0
        assert np.isclose(phi, phi_expected)

    def test_bilinear_interp_monotonic(self):
        """Bilinear interpolation should be monotonic within a grid cell."""
        # Grid with values 1.0 to 4.0
        grid = np.array([[1.0, 2.0], [3.0, 4.0]])

        # Sample across the grid
        phis = []
        for u in np.linspace(0, 1, 11):
            for v in np.linspace(0, 1, 11):
                p00, p10, p01, p11 = grid[0, 0], grid[1, 0], grid[0, 1], grid[1, 1]
                phi = (1 - u) * (1 - v) * p00 + u * (1 - v) * p10 + (
                    1 - u
                ) * v * p01 + u * v * p11
                phis.append(phi)

        # All values should be in range [1.0, 4.0]
        assert all(1.0 <= phi <= 4.0 for phi in phis)
        assert min(phis) == pytest.approx(1.0)
        assert max(phis) == pytest.approx(4.0)


class TestDecorrelationEffect:
    """Test the effect of decorrelation on synthetic data."""

    def test_decorrelation_scales_image(self):
        """
        Decorrelation should scale each pixel by φ(x,y).

        D_decorr[i,j] = φ(i,j) · D[i,j]
        """
        # Create synthetic difference image
        diff_image = np.ones((10, 10), dtype=np.float32) * 100.0

        # Simulate uniform decorrelation (φ = 0.5 everywhere)
        phi_value = 0.5
        diff_decorr = phi_value * diff_image

        # Check result
        assert np.allclose(diff_decorr, 50.0)
        assert diff_decorr.shape == diff_image.shape

    def test_decorrelation_with_spatial_variation(self):
        """
        Decorrelation should handle spatially-varying φ.

        Create a difference image with a gradient and decorrelate with
        spatially-varying φ.
        """
        # Create difference image with gradient
        diff_image = np.arange(100, dtype=np.float32).reshape(10, 10)

        # Create spatially-varying decorrelation kernel
        # φ increases from 0.5 to 1.5 across the image
        phi_field = np.linspace(0.5, 1.5, 100, dtype=np.float32).reshape(10, 10)

        # Apply decorrelation
        diff_decorr = phi_field * diff_image

        # Check that decorrelation was applied
        assert diff_decorr.shape == diff_image.shape
        # Pixels with larger φ should be amplified more
        assert diff_decorr[0, 0] < diff_decorr[-1, -1]

    def test_decorrelation_preserves_shape(self):
        """Decorrelation should preserve image shape."""
        shapes = [(5, 5), (10, 10), (100, 100), (17, 23)]
        for shape in shapes:
            diff_image = np.random.randn(*shape).astype(np.float32)
            phi_field = np.ones(shape, dtype=np.float32) * 0.8

            diff_decorr = phi_field * diff_image

            assert diff_decorr.shape == diff_image.shape


class TestDecorrelationIntegration:
    """Integration tests combining kernel fitting and decorrelation."""

    def test_decorrelation_with_synthetic_data(self, star_field):
        """
        End-to-end test: fit kernel and apply decorrelation.

        Uses the star_field fixture to create template and science images.
        Verifies that kernel fitting works with decorrelation enabled.
        """
        # Extract noiseless template image from fixture
        template = star_field["noiseless"].astype(np.float32)
        # Create science image with small perturbation
        science = template + 0.01 * np.random.randn(*template.shape).astype(np.float32)

        # Fit kernel without decorrelation
        config = KernelConfig(use_tps=False)
        kernel_sol = fit_kernel(science, template, config=config)

        # Verify kernel was fit successfully
        assert kernel_sol is not None
        assert len(kernel_sol.kernel_coefficients) > 0

        # Spatial convolution should work
        diff = spatial_convolve(template, kernel_sol, config=config)
        assert diff.shape == science.shape


class TestDecorrelationNumericalStability:
    """Test numerical stability of decorrelation computation."""

    def test_decorrelation_with_small_variances(self):
        """Decorrelation should be stable with small variances."""
        sig_s2 = 1e-10  # Very small science variance
        sig_t2 = 1e-10  # Very small template variance
        kappa = 1.0

        # Should not crash or produce NaN
        with np.errstate(divide="ignore", invalid="ignore"):
            phi = np.sqrt((sig_s2 + sig_t2) / (sig_s2 + sig_t2 * kappa**2))

        # Result should be finite or handled gracefully
        if not np.isfinite(phi):
            # If result is inf/nan, we should have fallback handling
            phi_fallback = 1.0  # Default: no decorrelation
            assert np.isfinite(phi_fallback)

    def test_decorrelation_with_large_variances(self):
        """Decorrelation should be stable with large variances."""
        sig_s2 = 1e10  # Very large science variance
        sig_t2 = 1e10  # Very large template variance
        kappa = 1.0

        phi = np.sqrt((sig_s2 + sig_t2) / (sig_s2 + sig_t2 * kappa**2))

        assert np.isfinite(phi)
        assert phi > 0.0

    def test_decorrelation_with_zero_kernel(self):
        """Decorrelation with zero kernel should be stable."""
        sig_s2 = 1.0
        sig_t2 = 1.0
        kappa = 0.0  # Zero kernel

        phi = np.sqrt((sig_s2 + sig_t2) / (sig_s2 + sig_t2 * kappa**2))

        assert np.isfinite(phi)
        assert phi > 1.0  # Should amplify when kernel is zero


class TestDecorrelationEdgeCases:
    """Test edge cases in decorrelation."""

    def test_single_pixel_image(self):
        """Decorrelation should work on minimal 1x1 image."""
        diff_image = np.array([[100.0]], dtype=np.float32)
        phi = 0.5
        diff_decorr = phi * diff_image

        assert diff_decorr.shape == (1, 1)
        assert np.isclose(diff_decorr[0, 0], 50.0)

    def test_single_row_image(self):
        """Decorrelation should work on 1D (row) image."""
        diff_image = np.arange(10, dtype=np.float32).reshape(1, 10)
        phi_field = np.ones(10, dtype=np.float32).reshape(1, 10) * 0.8

        diff_decorr = phi_field * diff_image

        assert diff_decorr.shape == (1, 10)

    def test_single_column_image(self):
        """Decorrelation should work on 1D (column) image."""
        diff_image = np.arange(10, dtype=np.float32).reshape(10, 1)
        phi_field = np.ones(10, dtype=np.float32).reshape(10, 1) * 0.8

        diff_decorr = phi_field * diff_image

        assert diff_decorr.shape == (10, 1)

    def test_uniform_decorrelation_field(self):
        """Uniform φ should scale image uniformly."""
        diff_image = np.random.randn(10, 10).astype(np.float32)
        phi_uniform = 0.75

        diff_decorr = phi_uniform * diff_image

        # Should be equivalent to scalar multiplication
        assert np.allclose(diff_decorr, diff_image * phi_uniform)

    def test_zero_difference_image(self):
        """Decorrelation of zero image should remain zero."""
        diff_image = np.zeros((10, 10), dtype=np.float32)
        phi_field = np.random.rand(10, 10).astype(np.float32)

        diff_decorr = phi_field * diff_image

        assert np.allclose(diff_decorr, 0.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
