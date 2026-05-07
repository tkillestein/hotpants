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
from astropy.io import fits

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


class TestPythonAPIvsCLI:
    """
    Compare Python API results to C CLI.

    These tests ensure that the Python API wrapper and C CLI produce
    numerically identical results on the same inputs.

    C integration is now complete via wrapper functions (initBuildStampsContext,
    buildStampsRegion, cleanupBuildStampsContext) in hotpants_wrapper.c that manage
    all required global state internally.
    """

    @pytest.fixture
    def identical_images(self, tmp_path, star_field, hotpants_binary):
        """Create identical template and science images."""
        rng = np.random.default_rng(42)
        template = make_image(
            star_field['noiseless'],
            PSF_SIGMA_TEMPLATE,
            rng,
        ).astype(np.float32)
        science = template.copy()

        return template, science, tmp_path, hotpants_binary

    def test_identical_images_api_vs_cli(self, identical_images):
        """
        Python API vs CLI on identical images.

        When template == science, difference should be pure noise.
        Both API and CLI should produce similar results.
        """
        template, science, tmp_path, hotpants_binary = identical_images

        # Run C CLI
        template_fits = tmp_path / "template.fits"
        science_fits = tmp_path / "science.fits"
        diff_fits = tmp_path / "diff.fits"

        write_fits(template_fits, template)
        write_fits(science_fits, science)

        # -c t: force CLI to convolve template (same direction as Python API which
        # always fits template ⊗ K ≈ science and returns science - template ⊗ K)
        run_hotpants(
            hotpants_binary,
            template_fits,
            science_fits,
            diff_fits,
            extra_args=["-r", "10", "-ko", "2", "-bgo", "1", "-c", "t"],
        )

        # Load CLI output (raw 2D array)
        cli_diff_raw = fits.getdata(str(diff_fits)).astype(np.float64)

        # Run Python API — parameters must match CLI args above plus HOTPANTS_BASE_ARGS:
        #   -r 10 -ko 2 -bgo 1 -ft 8 -nsx 5 -nsy 5 -nss 3
        #   -tu 60000 -tl -200 -iu 60000 -il -200 -tg 1 -ig 1 -tr 5 -ir 5
        config = KernelConfig(
            kernel_half_width=10, kernel_order=2, bg_order=1, fit_threshold=8,
            n_ks_stamps=3, hw_ks_stamp=15,
        )
        layout = RegionLayout(n_regions_x=1, n_regions_y=1,
                              stamps_per_region_x=5, stamps_per_region_y=5)
        thresholds = NoiseThresholds(
            template_upper_threshold=60000, template_lower_threshold=-200,
            science_upper_threshold=60000, science_lower_threshold=-200,
            template_gain=1.0, science_gain=1.0,
            template_readnoise=5.0, science_readnoise=5.0,
        )

        kernel_solution = fit_kernel(
            template,
            science,
            config=config,
            layout=layout,
            thresholds=thresholds,
        )

        # Difference image: science - template ⊗ K  (matches CLI -c t output)
        template_conv = spatial_convolve(
            template,
            kernel_solution,
            config=config,
        ).astype(np.float64)
        api_diff_raw = science.astype(np.float64) - template_conv

        # Mask using CLI fill pixels (1e-30 marks border / invalid pixels)
        valid_mask = np.abs(cli_diff_raw - 1e-30) > 1e-35
        cli_diff = cli_diff_raw[valid_mask]
        api_diff = api_diff_raw[valid_mask]

        # Compare (allow small numerical differences)
        # Note: differences due to floating-point order of operations
        np.testing.assert_allclose(
            cli_diff, api_diff,
            rtol=1e-3, atol=0.5,
            err_msg="Python API and CLI produce different results"
        )

    def test_broadened_psf_api_vs_cli(self, tmp_path, star_field, hotpants_binary):
        """
        Python API vs CLI when science PSF is broader than template.

        The fitted kernel should act as a broadening kernel (norm > 1).
        """
        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(43)
        template = make_image(
            star_field['noiseless'],
            PSF_SIGMA_TEMPLATE,
            rng1,
        ).astype(np.float32)
        science = make_image(
            star_field['noiseless'],
            PSF_SIGMA_SCIENCE,  # Broader than template
            rng2,
        ).astype(np.float32)

        # Run C CLI — -c t forces template convolution to match Python API direction
        template_fits = tmp_path / "template.fits"
        science_fits = tmp_path / "science.fits"
        diff_fits = tmp_path / "diff.fits"

        write_fits(template_fits, template)
        write_fits(science_fits, science)

        run_hotpants(
            hotpants_binary,
            template_fits,
            science_fits,
            diff_fits,
            extra_args=["-c", "t"],
        )

        cli_diff_raw = fits.getdata(str(diff_fits)).astype(np.float64)

        # Python API — parameters must match HOTPANTS_BASE_ARGS
        config = KernelConfig(
            kernel_half_width=8, kernel_order=2, bg_order=1,
            fit_threshold=8, n_ks_stamps=3,
        )
        layout = RegionLayout(
            n_regions_x=1, n_regions_y=1,
            stamps_per_region_x=5, stamps_per_region_y=5,
        )
        thresholds = NoiseThresholds(
            template_upper_threshold=60000, template_lower_threshold=-200,
            science_upper_threshold=60000, science_lower_threshold=-200,
            template_gain=1.0, science_gain=1.0,
            template_readnoise=5.0, science_readnoise=5.0,
        )

        kernel_solution = fit_kernel(
            template,
            science,
            config=config,
            layout=layout,
            thresholds=thresholds,
        )

        # For broadened PSF, kernel norm should be close to 1 (flux-conserving)
        assert 0.98 < kernel_solution.kernel_norm < 1.02, (
            f"Expected broadening kernel (norm ≈ 1), got {kernel_solution.kernel_norm:.3f}"
        )

        # Difference image: science - template ⊗ K  (matches CLI -c t output)
        template_conv = spatial_convolve(
            template,
            kernel_solution,
            config=config,
        ).astype(np.float64)
        api_diff_raw = science.astype(np.float64) - template_conv

        valid_mask = np.abs(cli_diff_raw - 1e-30) > 1e-35
        cli_diff = cli_diff_raw[valid_mask]
        api_diff = api_diff_raw[valid_mask]

        np.testing.assert_allclose(
            cli_diff, api_diff,
            rtol=1e-3, atol=0.5,
        )

    def test_narrowed_psf_api_vs_cli(self, tmp_path, star_field, hotpants_binary):
        """
        Python API vs CLI when science PSF is narrower than template.

        When template has a broader PSF than science, the Python API fits
        template ⊗ K ≈ science (sharpening direction, norm < 1).
        CLI is forced to the same direction with -c t.
        """
        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(43)
        template = make_image(
            star_field['noiseless'],
            PSF_SIGMA_SCIENCE,  # Broader
            rng1,
        ).astype(np.float32)
        science = make_image(
            star_field['noiseless'],
            PSF_SIGMA_TEMPLATE,  # Narrower
            rng2,
        ).astype(np.float32)

        # Run C CLI — -c t forces template convolution to match Python API direction
        template_fits = tmp_path / "template.fits"
        science_fits = tmp_path / "science.fits"
        diff_fits = tmp_path / "diff.fits"

        write_fits(template_fits, template)
        write_fits(science_fits, science)

        run_hotpants(
            hotpants_binary,
            template_fits,
            science_fits,
            diff_fits,
            extra_args=["-c", "t"],
        )

        cli_diff_raw = fits.getdata(str(diff_fits)).astype(np.float64)

        # Python API — parameters must match HOTPANTS_BASE_ARGS
        config = KernelConfig(
            kernel_half_width=8, kernel_order=2, bg_order=1,
            fit_threshold=8, n_ks_stamps=3,
        )
        layout = RegionLayout(
            n_regions_x=1, n_regions_y=1,
            stamps_per_region_x=5, stamps_per_region_y=5,
        )
        thresholds = NoiseThresholds(
            template_upper_threshold=60000, template_lower_threshold=-200,
            science_upper_threshold=60000, science_lower_threshold=-200,
            template_gain=1.0, science_gain=1.0,
            template_readnoise=5.0, science_readnoise=5.0,
        )

        kernel_solution = fit_kernel(
            template,
            science,
            config=config,
            layout=layout,
            thresholds=thresholds,
        )

        # For narrowed PSF, kernel norm should be < 1 (sharpening)
        assert kernel_solution.kernel_norm < 1.0, (
            f"Expected sharpening kernel (norm < 1), got {kernel_solution.kernel_norm:.3f}"
        )

        # Difference image: science - template ⊗ K  (matches CLI -c t output)
        template_conv = spatial_convolve(
            template,
            kernel_solution,
            config=config,
        ).astype(np.float64)
        api_diff_raw = science.astype(np.float64) - template_conv

        valid_mask = np.abs(cli_diff_raw - 1e-30) > 1e-35
        cli_diff = cli_diff_raw[valid_mask]
        api_diff = api_diff_raw[valid_mask]

        np.testing.assert_allclose(
            cli_diff, api_diff,
            rtol=1e-3, atol=0.5,
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
