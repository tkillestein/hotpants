"""
Unit tests for HOTPANTS Python API.

Tests cover:
- Configuration dataclass validation
- Input image validation
- API function signatures and error handling
- Basic end-to-end workflow
"""

import pytest
import numpy as np
from hotpants import (
    fit_kernel,
    spatial_convolve,
    KernelConfig,
    RegionLayout,
    NoiseThresholds,
    KernelSolution,
)


# =====================================================================
# Configuration Validation Tests
# =====================================================================

class TestKernelConfig:
    """Test KernelConfig dataclass validation."""

    def test_default_config(self):
        """Default config should be valid."""
        config = KernelConfig()
        assert config.kernel_half_width == 15
        assert config.kernel_order == 2
        assert config.bg_order == 1
        assert config.hw_ks_stamp == 10  # Must be <= kernel_half_width

    def test_custom_config(self):
        """Custom config with valid parameters."""
        config = KernelConfig(kernel_half_width=15, kernel_order=3, bg_order=2)
        assert config.kernel_half_width == 15
        assert config.kernel_order == 3

    def test_negative_kernel_half_width(self):
        """Negative kernel_half_width should raise ValueError."""
        with pytest.raises(ValueError, match="kernel_half_width must be positive"):
            KernelConfig(kernel_half_width=-1)

    def test_hw_ks_stamp_too_large(self):
        """hw_ks_stamp > kernel_half_width should raise ValueError."""
        with pytest.raises(ValueError, match="hw_ks_stamp must be <="):
            KernelConfig(kernel_half_width=10, hw_ks_stamp=15)

    def test_invalid_scale_fit_threshold(self):
        """scale_fit_threshold outside (0, 1] should raise ValueError."""
        with pytest.raises(ValueError, match="scale_fit_threshold"):
            KernelConfig(scale_fit_threshold=1.5)
        with pytest.raises(ValueError, match="scale_fit_threshold"):
            KernelConfig(scale_fit_threshold=0.0)


class TestRegionLayout:
    """Test RegionLayout dataclass validation."""

    def test_default_layout(self):
        """Default layout should be valid."""
        layout = RegionLayout()
        assert layout.n_regions_x == 1
        assert layout.n_regions_y == 1
        assert layout.stamps_per_region_x == 10
        assert layout.stamps_per_region_y == 10

    def test_invalid_regions(self):
        """Invalid region counts should raise ValueError."""
        with pytest.raises(ValueError, match="n_regions_"):
            RegionLayout(n_regions_x=0)


class TestNoiseThresholds:
    """Test NoiseThresholds dataclass validation."""

    def test_default_thresholds(self):
        """Default thresholds should be valid."""
        thresh = NoiseThresholds()
        assert thresh.template_upper_threshold > thresh.template_lower_threshold

    def test_invalid_gain(self):
        """Non-positive gain should raise ValueError."""
        with pytest.raises(ValueError, match="gain"):
            NoiseThresholds(template_gain=0.0)

    def test_invalid_threshold_order(self):
        """upper_threshold <= lower_threshold should raise ValueError."""
        with pytest.raises(ValueError, match="upper_threshold"):
            NoiseThresholds(
                template_upper_threshold=100.0,
                template_lower_threshold=200.0
            )


# =====================================================================
# Input Validation Tests
# =====================================================================

class TestInputValidation:
    """Test input image validation."""

    @pytest.fixture
    def valid_images(self):
        """Create valid test images."""
        shape = (64, 64)
        template = np.random.randn(*shape).astype(np.float32)
        science = np.random.randn(*shape).astype(np.float32)
        return template, science

    def test_wrong_dtype_template(self, valid_images):
        """float64 template should raise TypeError."""
        template, science = valid_images
        template = template.astype(np.float64)
        with pytest.raises(TypeError, match="float32"):
            fit_kernel(template, science)

    def test_wrong_dtype_science(self, valid_images):
        """float64 science should raise TypeError."""
        template, science = valid_images
        science = science.astype(np.float64)
        with pytest.raises(TypeError, match="float32"):
            fit_kernel(template, science)

    def test_shape_mismatch(self, valid_images):
        """Mismatched image shapes should raise ValueError."""
        template, science = valid_images
        science = science[:-1, :-1]  # Reduce size
        with pytest.raises(ValueError, match="shape mismatch"):
            fit_kernel(template, science)

    def test_nan_in_template(self, valid_images):
        """NaN in template should raise ValueError."""
        template, science = valid_images
        template[0, 0] = np.nan
        with pytest.raises(ValueError, match="NaN or Inf"):
            fit_kernel(template, science)

    def test_inf_in_science(self, valid_images):
        """Inf in science should raise ValueError."""
        template, science = valid_images
        science[0, 0] = np.inf
        with pytest.raises(ValueError, match="NaN or Inf"):
            fit_kernel(template, science)

    def test_1d_array(self):
        """1D array should raise ValueError."""
        image_1d = np.random.randn(64).astype(np.float32)
        science_2d = np.random.randn(64, 64).astype(np.float32)
        with pytest.raises(ValueError, match="2D array"):
            fit_kernel(image_1d, science_2d)

    def test_optional_noise_validation(self, valid_images):
        """Noise image validation should work."""
        template, science = valid_images
        noise = np.random.randn(*template.shape).astype(np.float64)  # Wrong dtype
        with pytest.raises(TypeError, match="float32"):
            fit_kernel(template, science, noise=noise)

    def test_optional_noise_shape_mismatch(self, valid_images):
        """Noise shape mismatch should raise ValueError."""
        template, science = valid_images
        noise = np.random.randn(32, 32).astype(np.float32)  # Wrong shape
        with pytest.raises(ValueError, match="Noise shape"):
            fit_kernel(template, science, noise=noise)


# =====================================================================
# API Function Tests
# =====================================================================

class TestFitKernel:
    """Test fit_kernel function."""

    @pytest.fixture
    def valid_images(self):
        """Create valid test images."""
        np.random.seed(42)
        shape = (128, 128)
        template = np.random.randn(*shape).astype(np.float32)
        science = template + np.random.randn(*shape).astype(np.float32) * 0.1
        return template, science

    def test_basic_fit_kernel(self, valid_images):
        """Basic kernel fitting should work."""
        template, science = valid_images
        result = fit_kernel(template, science)

        assert isinstance(result, KernelSolution)
        assert result.kernel_coefficients is not None
        assert result.chi2 >= 0
        assert result.kernel_norm > 0

    def test_fit_kernel_with_config(self, valid_images):
        """fit_kernel with custom config should work."""
        template, science = valid_images
        config = KernelConfig(kernel_half_width=15, kernel_order=3)
        result = fit_kernel(template, science, config=config)

        assert isinstance(result, KernelSolution)

    def test_fit_kernel_with_layout(self, valid_images):
        """fit_kernel with custom layout should work."""
        template, science = valid_images
        layout = RegionLayout(n_regions_x=1, n_regions_y=1, stamps_per_region_x=5)
        result = fit_kernel(template, science, layout=layout)

        assert isinstance(result, KernelSolution)

    def test_fit_kernel_with_thresholds(self, valid_images):
        """fit_kernel with custom thresholds should work."""
        template, science = valid_images
        thresholds = NoiseThresholds(
            template_gain=2.0,
            science_gain=2.0,
            template_readnoise=10.0
        )
        result = fit_kernel(template, science, thresholds=thresholds)

        assert isinstance(result, KernelSolution)

    def test_fit_kernel_with_noise(self, valid_images):
        """fit_kernel with noise image should work."""
        template, science = valid_images
        noise = np.ones_like(template) * 0.1
        result = fit_kernel(template, science, noise=noise)

        assert isinstance(result, KernelSolution)

    def test_fit_kernel_verbose(self, valid_images, capsys):
        """fit_kernel with verbose flag should not crash."""
        template, science = valid_images
        result = fit_kernel(template, science, verbose=1)

        assert isinstance(result, KernelSolution)

    def test_kernel_solution_repr(self, valid_images):
        """KernelSolution repr should be readable."""
        template, science = valid_images
        result = fit_kernel(template, science)
        repr_str = repr(result)

        assert "KernelSolution" in repr_str
        assert "chi2" in repr_str


class TestSpatialConvolve:
    """Test spatial_convolve function."""

    @pytest.fixture
    def kernel_solution(self):
        """Create a dummy kernel solution."""
        return KernelSolution(
            chi2=0.0,
            kernel_norm=1.0,
            mean_sigma=10.0,
            scatter_sigma=1.0,
            n_skipped_stamps=0,
            kernel_coefficients=np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        )

    def test_basic_convolve(self, kernel_solution):
        """Basic convolution should work."""
        image = np.random.randn(64, 64).astype(np.float32)
        output = spatial_convolve(image, kernel_solution)

        assert isinstance(output, np.ndarray)
        assert output.shape == image.shape
        assert output.dtype == np.float32

    def test_convolve_with_config(self, kernel_solution):
        """Convolution with custom config should work."""
        image = np.random.randn(64, 64).astype(np.float32)
        config = KernelConfig(kernel_half_width=15)
        output = spatial_convolve(image, kernel_solution, config=config)

        assert output.shape == image.shape

    def test_convolve_with_preallocated_output(self, kernel_solution):
        """Convolution with pre-allocated output should work."""
        image = np.random.randn(64, 64).astype(np.float32)
        output = np.zeros_like(image)
        result = spatial_convolve(image, kernel_solution, output=output)

        assert result is output
        assert output.shape == image.shape

    def test_convolve_output_shape_mismatch(self, kernel_solution):
        """Pre-allocated output with wrong shape should raise ValueError."""
        image = np.random.randn(64, 64).astype(np.float32)
        output = np.zeros((32, 32), dtype=np.float32)

        with pytest.raises(ValueError, match="shape"):
            spatial_convolve(image, kernel_solution, output=output)

    def test_convolve_wrong_dtype(self, kernel_solution):
        """Wrong dtype image should raise TypeError."""
        image = np.random.randn(64, 64).astype(np.float64)

        with pytest.raises(TypeError, match="float32"):
            spatial_convolve(image, kernel_solution)


# =====================================================================
# Integration Tests
# =====================================================================

class TestIntegration:
    """Integration tests for full workflow."""

    def test_fit_and_convolve_workflow(self):
        """Full workflow: fit kernel and apply to image."""
        np.random.seed(42)
        shape = (128, 128)

        # Create synthetic images
        template = np.random.randn(*shape).astype(np.float32) * 100 + 1000
        science = template + np.random.randn(*shape).astype(np.float32) * 10

        # Fit kernel
        config = KernelConfig()
        layout = RegionLayout()
        thresholds = NoiseThresholds()

        kernel_solution = fit_kernel(
            template, science,
            config=config,
            layout=layout,
            thresholds=thresholds
        )

        # Apply kernel
        difference = spatial_convolve(science, kernel_solution, config=config)

        # Verify outputs
        assert difference.shape == science.shape
        assert difference.dtype == np.float32
        assert kernel_solution.kernel_norm > 0

    def test_identical_images_fit(self):
        """Fitting kernel on identical images should give near-identity."""
        np.random.seed(42)
        shape = (96, 96)
        template = np.random.randn(*shape).astype(np.float32)
        science = template.copy()

        result = fit_kernel(template, science)

        # For identical images, kernel should be near delta (norm ≈ 1)
        # and chi2 should be small
        assert 0.5 < result.kernel_norm < 2.0  # Some tolerance

    def test_c_contiguous_conversion(self):
        """Non-C-contiguous arrays should be converted."""
        np.random.seed(42)
        shape = (64, 64)

        # Create Fortran-contiguous (column-major) arrays
        template = np.asfortranarray(np.random.randn(*shape).astype(np.float32))
        science = np.asfortranarray(np.random.randn(*shape).astype(np.float32))

        assert not template.flags['C_CONTIGUOUS']
        assert not science.flags['C_CONTIGUOUS']

        # Should work despite non-C-contiguous input
        result = fit_kernel(template, science)
        assert isinstance(result, KernelSolution)
