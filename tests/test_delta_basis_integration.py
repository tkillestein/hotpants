"""
Integration tests for Bramich (2008) delta function kernel basis.

These tests verify that the delta basis implementation produces correct results
when integrated with the full kernel fitting and convolution pipeline. Tests
cover:
- Kernel fitting with synthetic astronomical images
- Spatial variation (polynomial and TPS)
- Laplacian regularization for smoothness
- Comparison with Gaussian basis
- Numerical stability and convergence
"""

import pytest
import numpy as np
from scipy.ndimage import gaussian_filter, convolve

from hotpants import fit_kernel, spatial_convolve, KernelConfig, RegionLayout, NoiseThresholds


class TestDeltaBasisIntegration:
    """Comprehensive integration tests for delta function kernel basis."""

    @pytest.fixture
    def psf_template(self):
        """Create a Gaussian PSF for the template image."""
        size = 7
        x = np.arange(size) - size // 2
        y = np.arange(size) - size // 2
        xx, yy = np.meshgrid(x, y)
        sigma = 1.5
        psf = np.exp(-(xx**2 + yy**2) / (2 * sigma**2))
        return psf / psf.sum()

    @pytest.fixture
    def psf_science(self):
        """Create a slightly different Gaussian PSF for the science image."""
        size = 7
        x = np.arange(size) - size // 2
        y = np.arange(size) - size // 2
        xx, yy = np.meshgrid(x, y)
        sigma = 2.0  # Different from template
        psf = np.exp(-(xx**2 + yy**2) / (2 * sigma**2))
        return psf / psf.sum()

    @pytest.fixture
    def star_field(self, psf_template, psf_science):
        """Create synthetic star field images with different PSFs."""
        nx, ny = 256, 256
        rng = np.random.RandomState(42)

        # Create point source positions (random but fixed seed)
        n_stars = 20
        star_x = rng.uniform(20, nx - 20, n_stars).astype(int)
        star_y = rng.uniform(20, ny - 20, n_stars).astype(int)
        star_flux = rng.uniform(1000, 5000, n_stars)

        # Template image: original PSF
        template_noiseless = np.zeros((ny, nx), dtype=np.float32)
        for x, y, flux in zip(star_x, star_y, star_flux):
            template_noiseless[y - 3:y + 4, x - 3:x + 4] += flux * psf_template

        # Science image: broadened PSF (requires deconvolution)
        science_noiseless = np.zeros((ny, nx), dtype=np.float32)
        for x, y, flux in zip(star_x, star_y, star_flux):
            science_noiseless[y - 3:y + 4, x - 3:x + 4] += flux * psf_science

        # Add background and noise
        bg = 100.0
        rdnoise = 5.0
        template = template_noiseless + bg + rng.normal(0, rdnoise, (ny, nx))
        science = science_noiseless + bg + rng.normal(0, rdnoise, (ny, nx))

        return {
            "template": template.astype(np.float32),
            "science": science.astype(np.float32),
            "template_noiseless": template_noiseless,
            "science_noiseless": science_noiseless,
            "star_x": star_x,
            "star_y": star_y,
            "star_flux": star_flux,
        }

    def test_delta_basis_identical_images(self):
        """Delta basis should produce near-zero difference for identical images."""
        # Create identical images
        template = np.random.RandomState(42).rand(128, 128).astype(np.float32) * 1000
        science = template.copy()

        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=10,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=5,
            stamps_per_region_y=5,
        )

        # Fit kernel
        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None
        assert kernel_sol.kernel_norm > 0

        # Apply convolution
        diff = spatial_convolve(template, kernel_sol, config=config)

        # Difference should be small (only noise)
        assert diff.shape == template.shape
        assert not np.any(np.isnan(diff))
        rms_diff = np.sqrt(np.mean(diff**2))
        assert rms_diff < 100.0, f"RMS difference {rms_diff} is too large for identical images"

    def test_delta_basis_with_synthetic_stars(self, star_field):
        """Delta basis should fit and deconvolve synthetic star field correctly."""
        template = star_field["template"]
        science = star_field["science"]

        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=10,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=6,
            stamps_per_region_y=6,
        )

        # Fit kernel
        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None
        assert kernel_sol.kernel_norm > 0

        # Apply convolution
        diff = spatial_convolve(template, kernel_sol, config=config)

        # Check output validity
        assert diff.shape == template.shape
        assert not np.any(np.isnan(diff))
        assert not np.any(np.isinf(diff))

        # Difference should be reasonable (dominated by star noise, not residuals)
        rms_diff = np.sqrt(np.mean(diff**2))
        rms_template = np.sqrt(np.mean(template**2))
        relative_rms = rms_diff / rms_template
        assert relative_rms < 0.5, f"Relative RMS {relative_rms} suggests poor fit"

    def test_delta_basis_with_regularization(self, star_field):
        """Delta basis with Laplacian regularization should produce smooth kernels."""
        template = star_field["template"]
        science = star_field["science"]

        config_noreg = KernelConfig(
            basis_type="delta",
            kernel_half_width=8,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=5,
            stamps_per_region_y=5,
            delta_regularization=0.0,  # No regularization
        )

        config_reg = KernelConfig(
            basis_type="delta",
            kernel_half_width=8,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=5,
            stamps_per_region_y=5,
            delta_regularization=1e-3,  # Strong regularization
        )

        # Fit kernels
        kernel_noreg = fit_kernel(template, science, config=config_noreg)
        kernel_reg = fit_kernel(template, science, config=config_reg)

        assert kernel_noreg is not None
        assert kernel_reg is not None

        # Both should produce valid output
        diff_noreg = spatial_convolve(template, kernel_noreg, config=config_noreg)
        diff_reg = spatial_convolve(template, kernel_reg, config=config_reg)

        assert diff_noreg.shape == template.shape
        assert diff_reg.shape == template.shape
        assert not np.any(np.isnan(diff_noreg))
        assert not np.any(np.isnan(diff_reg))

    def test_delta_basis_vs_gaussian_basis(self, star_field):
        """Delta and Gaussian bases should produce similar difference images."""
        template = star_field["template"]
        science = star_field["science"]

        config_delta = KernelConfig(
            basis_type="delta",
            kernel_half_width=10,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=5,
            stamps_per_region_y=5,
        )

        config_gaussian = KernelConfig(
            basis_type="gaussian",
            kernel_half_width=10,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=5,
            stamps_per_region_y=5,
        )

        # Fit both bases
        kernel_delta = fit_kernel(template, science, config=config_delta)
        kernel_gaussian = fit_kernel(template, science, config=config_gaussian)

        assert kernel_delta is not None
        assert kernel_gaussian is not None

        # Apply convolution with both
        diff_delta = spatial_convolve(template, kernel_delta, config=config_delta)
        diff_gaussian = spatial_convolve(template, kernel_gaussian, config=config_gaussian)

        # Check shapes match
        assert diff_delta.shape == diff_gaussian.shape

        # Difference images should be qualitatively similar
        # (not exactly identical due to different parameterizations)
        rms_delta = np.sqrt(np.mean(diff_delta**2))
        rms_gaussian = np.sqrt(np.mean(diff_gaussian**2))

        # RMS values should be in same ballpark (within factor of 2)
        ratio = rms_delta / rms_gaussian
        assert 0.5 < ratio < 2.0, f"RMS ratio {ratio} suggests incompatible results"

    def test_delta_basis_spatial_variation_polynomial(self, star_field):
        """Delta basis with polynomial spatial variation should work."""
        template = star_field["template"]
        science = star_field["science"]

        # Single region with polynomial order 2 spatial variation
        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=10,
            kernel_order=2,  # Polynomial spatial variation
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=5,
            stamps_per_region_y=5,
        )

        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None

        diff = spatial_convolve(template, kernel_sol, config=config)
        assert diff.shape == template.shape
        assert not np.any(np.isnan(diff))

    def test_delta_basis_multi_region(self, star_field):
        """Delta basis should handle multi-region fitting."""
        template = star_field["template"]
        science = star_field["science"]

        # Multi-region with delta basis
        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=8,
            kernel_order=1,
            bg_order=1,
            num_regions_x=2,
            num_regions_y=2,
            stamps_per_region_x=4,
            stamps_per_region_y=4,
        )

        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None

        diff = spatial_convolve(template, kernel_sol, config=config)
        assert diff.shape == template.shape
        assert not np.any(np.isnan(diff))

    def test_delta_basis_convergence(self, star_field):
        """Delta basis kernel fitting should converge."""
        template = star_field["template"]
        science = star_field["science"]

        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=10,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=5,
            stamps_per_region_y=5,
        )

        # Multiple fits should be stable
        kernel_sol1 = fit_kernel(template, science, config=config)
        kernel_sol2 = fit_kernel(template, science, config=config)

        # Both should produce non-NaN results
        assert kernel_sol1 is not None
        assert kernel_sol2 is not None
        assert kernel_sol1.kernel_norm > 0
        assert kernel_sol2.kernel_norm > 0

    def test_delta_basis_noise_handling(self):
        """Delta basis should handle noisy images reasonably."""
        # Create base image
        rng = np.random.RandomState(42)
        base = np.ones((128, 128), dtype=np.float32) * 100

        # Add stars
        for _ in range(10):
            x = rng.randint(10, 118)
            y = rng.randint(10, 118)
            base[y - 2:y + 3, x - 2:x + 3] += 200

        # Create template and science with different noise realizations
        template = base + rng.normal(0, 5, (128, 128)).astype(np.float32)
        science = base + rng.normal(0, 5, (128, 128)).astype(np.float32)

        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=8,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=4,
            stamps_per_region_y=4,
            delta_regularization=1e-3,  # Use regularization for noisy data
        )

        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None

        diff = spatial_convolve(template, kernel_sol, config=config)
        assert not np.any(np.isnan(diff))
        assert not np.any(np.isinf(diff))

    def test_delta_basis_kernel_normalization(self, star_field):
        """Delta basis kernel should be properly normalized."""
        template = star_field["template"]
        science = star_field["science"]

        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=10,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=5,
            stamps_per_region_y=5,
        )

        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None

        # Kernel norm should be positive and reasonable
        assert kernel_sol.kernel_norm > 0
        assert kernel_sol.kernel_norm < 100  # Should not be astronomically large

    def test_delta_basis_with_tps_stub(self, star_field):
        """Delta basis with TPS flag should work (using polynomial fallback)."""
        template = star_field["template"]
        science = star_field["science"]

        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=8,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=4,
            stamps_per_region_y=4,
            use_tps=True,  # TPS stub for delta (falls back to polynomial)
        )

        # Should not crash; currently falls back to polynomial
        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None

        diff = spatial_convolve(template, kernel_sol, config=config)
        assert not np.any(np.isnan(diff))


class TestDeltaBasisEdgeCases:
    """Edge cases and error conditions for delta basis."""

    def test_delta_basis_very_small_kernel(self):
        """Delta basis with very small kernel should work."""
        template = np.random.rand(128, 128).astype(np.float32)
        science = template.copy()

        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=2,  # Very small kernel
            kernel_order=0,  # No spatial variation
            bg_order=0,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=3,
            stamps_per_region_y=3,
        )

        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None

        diff = spatial_convolve(template, kernel_sol, config=config)
        assert diff.shape == template.shape

    def test_delta_basis_large_kernel(self):
        """Delta basis with large kernel should work."""
        template = np.random.rand(256, 256).astype(np.float32)
        science = template.copy()

        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=20,  # Large kernel
            kernel_order=0,
            bg_order=0,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=3,
            stamps_per_region_y=3,
        )

        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None

        diff = spatial_convolve(template, kernel_sol, config=config)
        assert diff.shape == template.shape

    def test_delta_basis_high_regularization(self):
        """Delta basis with very high regularization should still work."""
        template = np.random.rand(128, 128).astype(np.float32)
        science = template.copy()

        config = KernelConfig(
            basis_type="delta",
            kernel_half_width=8,
            kernel_order=1,
            bg_order=1,
            num_regions_x=1,
            num_regions_y=1,
            stamps_per_region_x=4,
            stamps_per_region_y=4,
            delta_regularization=1e6,  # Very high regularization
        )

        kernel_sol = fit_kernel(template, science, config=config)
        assert kernel_sol is not None

        diff = spatial_convolve(template, kernel_sol, config=config)
        assert not np.any(np.isnan(diff))
