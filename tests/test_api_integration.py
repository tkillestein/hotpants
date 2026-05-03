"""
Integration tests: Python API vs C CLI comparison.

These tests verify that the Python API produces results consistent with
the C command-line interface. They use the same synthetic image fixtures
as the C regression suite.

Note: These tests will be skipped initially (until cffi C integration is
complete), but they document the expected behavior for the Python API.
"""

import pytest
import numpy as np

from .conftest import (
    NX, NY,
    make_image,
    run_hotpants,
    write_fits,
    load_diff,
    PSF_SIGMA_TEMPLATE,
    PSF_SIGMA_SCIENCE,
    PSF_SIGMA_EXTRA,
    BACKGROUND,
    RDNOISE,
)

# Try to import Python API; skip tests if not available
pytest.importorskip("hotpants")
from hotpants import fit_kernel, spatial_convolve, KernelConfig, RegionLayout, NoiseThresholds


@pytest.mark.skip(reason="C binding integration not yet complete")
class TestPythonAPIvsCLI:
    """
    Compare Python API results to C CLI.

    These tests ensure that the Python API wrapper and C CLI produce
    numerically identical results on the same inputs.
    """

    @pytest.fixture
    def identical_images(self, tmp_path):
        """Create identical template and science images."""
        template = make_image(
            seed=42,
            psf_sigma=PSF_SIGMA_TEMPLATE,
            background=BACKGROUND,
            rdnoise=RDNOISE,
        )
        science = template.copy()

        return template, science, tmp_path

    def test_identical_images_api_vs_cli(self, identical_images):
        """
        Python API vs CLI on identical images.

        When template == science, difference should be pure noise.
        Both API and CLI should produce similar results.
        """
        template, science, tmp_path = identical_images

        # Run C CLI
        template_fits = tmp_path / "template.fits"
        science_fits = tmp_path / "science.fits"
        diff_fits = tmp_path / "diff.fits"

        write_fits(template_fits, template)
        write_fits(science_fits, science)

        run_hotpants(
            str(template_fits), str(science_fits), str(diff_fits),
            hwKernel=10, kerOrder=2, bgOrder=1,
            nRegX=1, nRegY=1,
        )

        cli_diff = load_diff(diff_fits)

        # Run Python API
        config = KernelConfig(kernel_half_width=10, kernel_order=2, bg_order=1)
        layout = RegionLayout(n_regions_x=1, n_regions_y=1)
        thresholds = NoiseThresholds()

        kernel_solution = fit_kernel(
            template.astype(np.float32),
            science.astype(np.float32),
            config=config,
            layout=layout,
            thresholds=thresholds,
        )

        api_diff = spatial_convolve(
            science.astype(np.float32),
            kernel_solution,
            config=config,
        )

        # Compare (allow small numerical differences)
        # Note: differences due to floating-point order of operations
        np.testing.assert_allclose(
            cli_diff, api_diff,
            rtol=1e-4, atol=1e-3,
            err_msg="Python API and CLI produce different results"
        )

    def test_broadened_psf_api_vs_cli(self, tmp_path):
        """
        Python API vs CLI when science PSF is broader than template.

        The fitted kernel should act as a broadening kernel.
        """
        template = make_image(
            seed=42,
            psf_sigma=PSF_SIGMA_TEMPLATE,
            background=BACKGROUND,
            rdnoise=RDNOISE,
        )
        science = make_image(
            seed=43,
            psf_sigma=PSF_SIGMA_SCIENCE,  # Broader than template
            background=BACKGROUND,
            rdnoise=RDNOISE,
        )

        # Run C CLI
        template_fits = tmp_path / "template.fits"
        science_fits = tmp_path / "science.fits"
        diff_fits = tmp_path / "diff.fits"

        write_fits(template_fits, template)
        write_fits(science_fits, science)

        run_hotpants(
            str(template_fits), str(science_fits), str(diff_fits),
            hwKernel=10, kerOrder=2, bgOrder=1,
        )

        cli_diff = load_diff(diff_fits)

        # Run Python API
        config = KernelConfig()
        layout = RegionLayout()
        thresholds = NoiseThresholds()

        kernel_solution = fit_kernel(
            template.astype(np.float32),
            science.astype(np.float32),
            config=config,
            layout=layout,
            thresholds=thresholds,
        )

        api_diff = spatial_convolve(
            science.astype(np.float32),
            kernel_solution,
            config=config,
        )

        # For broadened PSF, kernel norm should be > 1 (broadening)
        # Both should produce similar results
        np.testing.assert_allclose(
            cli_diff, api_diff,
            rtol=1e-3, atol=1e-2,
        )

    def test_narrowed_psf_api_vs_cli(self, tmp_path):
        """
        Python API vs CLI when science PSF is narrower than template.

        The fitted kernel should act as a sharpening kernel.
        """
        template = make_image(
            seed=42,
            psf_sigma=PSF_SIGMA_SCIENCE,  # Broader
            background=BACKGROUND,
            rdnoise=RDNOISE,
        )
        science = make_image(
            seed=43,
            psf_sigma=PSF_SIGMA_TEMPLATE,  # Narrower
            background=BACKGROUND,
            rdnoise=RDNOISE,
        )

        # Run C CLI
        template_fits = tmp_path / "template.fits"
        science_fits = tmp_path / "science.fits"
        diff_fits = tmp_path / "diff.fits"

        write_fits(template_fits, template)
        write_fits(science_fits, science)

        run_hotpants(
            str(template_fits), str(science_fits), str(diff_fits),
            hwKernel=10, kerOrder=2, bgOrder=1,
        )

        cli_diff = load_diff(diff_fits)

        # Run Python API
        config = KernelConfig()
        layout = RegionLayout()
        thresholds = NoiseThresholds()

        kernel_solution = fit_kernel(
            template.astype(np.float32),
            science.astype(np.float32),
            config=config,
            layout=layout,
            thresholds=thresholds,
        )

        api_diff = spatial_convolve(
            science.astype(np.float32),
            kernel_solution,
            config=config,
        )

        # For narrowed PSF, kernel norm should be < 1 (sharpening)
        np.testing.assert_allclose(
            cli_diff, api_diff,
            rtol=1e-3, atol=1e-2,
        )


class TestPythonAPIStructure:
    """
    Tests for Python API structure that don't require C integration.

    These tests verify the API layer works correctly without needing
    the C code to be fully implemented.
    """

    def test_kernel_config_to_globals(self):
        """KernelConfig should map to correct global variable names."""
        config = KernelConfig(
            kernel_half_width=15,
            kernel_order=3,
            bg_order=2,
        )

        # Verify config can be used to construct global dict
        global_dict = {
            'hwKernel': config.kernel_half_width,
            'kerOrder': config.kernel_order,
            'bgOrder': config.bg_order,
        }

        assert global_dict['hwKernel'] == 15
        assert global_dict['kerOrder'] == 3
        assert global_dict['bgOrder'] == 2

    def test_region_layout_to_globals(self):
        """RegionLayout should map to correct global variable names."""
        layout = RegionLayout(
            n_regions_x=2,
            n_regions_y=2,
            stamps_per_region_x=8,
        )

        global_dict = {
            'nRegX': layout.n_regions_x,
            'nRegY': layout.n_regions_y,
            'nStampX': layout.stamps_per_region_x,
            'nStampY': layout.stamps_per_region_y,
        }

        assert global_dict['nRegX'] == 2
        assert global_dict['nRegY'] == 2
        assert global_dict['nStampX'] == 8

    def test_noise_thresholds_to_globals(self):
        """NoiseThresholds should map to correct global variable names."""
        thresh = NoiseThresholds(
            template_upper_threshold=30000.0,
            template_gain=2.5,
            science_readnoise=15.0,
        )

        global_dict = {
            'tUThresh': thresh.template_upper_threshold,
            'tGain': thresh.template_gain,
            'iRdnoise': thresh.science_readnoise,
        }

        assert global_dict['tUThresh'] == 30000.0
        assert global_dict['tGain'] == 2.5
        assert global_dict['iRdnoise'] == 15.0
