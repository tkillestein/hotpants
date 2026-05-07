"""
Diagnostic test: Compare Python API vs C CLI for different PSF scenarios.

This script runs both API and CLI on identical test cases and logs intermediate
values to identify exactly where they diverge.
"""

import subprocess
from pathlib import Path
import tempfile

import numpy as np
import pytest
from astropy.io import fits

from .conftest import (
    NX, NY,
    make_image,
    run_hotpants,
    write_fits,
    PSF_SIGMA_TEMPLATE,
    PSF_SIGMA_SCIENCE,
    BACKGROUND,
    RDNOISE,
)

pytest.importorskip("hotpants")
from hotpants import fit_kernel, spatial_convolve, KernelConfig, RegionLayout, NoiseThresholds
from hotpants._hotpants_ffi import get_global_int


class TestDiagnosticPSF:
    """Diagnostic tests with detailed logging for different PSF scenarios."""

    @pytest.fixture
    def broadened_psf_images(self, star_field):
        """Create template (narrow PSF) and science (broad PSF) images."""
        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(43)
        template = make_image(
            star_field['noiseless'],
            PSF_SIGMA_TEMPLATE,
            rng1,
        ).astype(np.float32)
        science = make_image(
            star_field['noiseless'],
            PSF_SIGMA_SCIENCE,  # Broader
            rng2,
        ).astype(np.float32)
        return template, science

    def test_diagnostic_broadened_psf(self, broadened_psf_images, tmp_path, hotpants_binary, star_field):
        """
        Run broadened PSF test with detailed diagnostic logging.

        Outputs:
        - Intermediate values from both API and CLI
        - Side-by-side comparison of key metrics
        - CSV file with detailed results
        """
        template, science = broadened_psf_images

        # ============ Setup: Create test FITS files ============
        template_fits = tmp_path / "template.fits"
        science_fits = tmp_path / "science.fits"
        diff_fits_cli = tmp_path / "diff_cli.fits"
        diff_fits_api = tmp_path / "diff_api.fits"

        write_fits(template_fits, template)
        write_fits(science_fits, science)

        print("\n" + "="*80)
        print("DIAGNOSTIC TEST: Broadened PSF (Science σ=2.30 px, Template σ=1.5 px)")
        print("="*80)

        # ============ Run C CLI ============
        print("\n[1/3] Running C CLI...")
        run_hotpants(
            hotpants_binary,
            template_fits,
            science_fits,
            diff_fits_cli,
            extra_args=["-c", "t", "-v", "2"],  # -v 2 for verbose output
        )

        cli_diff_raw = fits.getdata(str(diff_fits_cli)).astype(np.float64)
        valid_mask_cli = np.abs(cli_diff_raw - 1e-30) > 1e-35
        cli_diff = cli_diff_raw[valid_mask_cli]

        print(f"  CLI diff image: shape={cli_diff_raw.shape}, valid pixels={valid_mask_cli.sum()}")
        print(f"  CLI diff stats: min={cli_diff.min():.2f}, max={cli_diff.max():.2f}, "
              f"mean={cli_diff.mean():.2f}, std={cli_diff.std():.2f}")

        # ============ Run Python API with logging ============
        print("\n[2/3] Running Python API...")

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

        # Fit kernel and capture intermediate values
        kernel_solution = fit_kernel(
            template,
            science,
            config=config,
            layout=layout,
            thresholds=thresholds,
        )

        print(f"  Kernel norm: {kernel_solution.kernel_norm:.6f}")
        print(f"  Mean sigma: {kernel_solution.mean_sigma:.6f}")
        print(f"  Scatter sigma: {kernel_solution.scatter_sigma:.6f}")
        print(f"  Chi2: {kernel_solution.chi2:.6f}")

        # Check nS (number of valid stamps)
        nS = get_global_int("nS")
        print(f"  Valid stamps (nS): {nS}")

        # Get kernel coefficients for comparison
        api_coeffs = kernel_solution.kernel_coefficients
        print(f"  Kernel coefficients shape: {api_coeffs.shape}")
        print(f"  First 5 coefficients: {api_coeffs[:5]}")

        # Compute difference image from API
        template_conv = spatial_convolve(
            template,
            kernel_solution,
            config=config,
        ).astype(np.float64)
        api_diff_raw = science.astype(np.float64) - template_conv

        valid_mask_api = np.abs(api_diff_raw - 1e-30) > 1e-35  # Use same mask logic
        api_diff = api_diff_raw[valid_mask_api]

        print(f"  API diff image: shape={api_diff_raw.shape}, valid pixels={valid_mask_api.sum()}")
        print(f"  API diff stats: min={api_diff.min():.2f}, max={api_diff.max():.2f}, "
              f"mean={api_diff.mean():.2f}, std={api_diff.std():.2f}")

        # ============ Compare outputs ============
        print("\n[3/3] Comparing outputs...")

        # Use the smaller valid set (should be same)
        min_valid = min(cli_diff.size, api_diff.size)
        if min_valid > 0:
            cli_subset = cli_diff[:min_valid]
            api_subset = api_diff[:min_valid]

            pixel_diffs = np.abs(cli_subset - api_subset)
            print(f"  Pixel diff stats: min={pixel_diffs.min():.6f}, max={pixel_diffs.max():.6f}, "
                  f"mean={pixel_diffs.mean():.6f}, std={pixel_diffs.std():.6f}")

            # Report if max diff exceeds tolerance
            if pixel_diffs.max() > 0.5:
                print(f"  ⚠️  MAX DIFF EXCEEDS TOLERANCE: {pixel_diffs.max():.2f} ADU > 0.5 ADU")

                # Find worst pixels
                worst_idx = np.argsort(pixel_diffs)[-5:]
                print(f"  Worst 5 pixels:")
                for idx in reversed(worst_idx):
                    print(f"    Pixel {idx}: CLI={cli_subset[idx]:.2f}, API={api_subset[idx]:.2f}, "
                          f"Diff={pixel_diffs[idx]:.2f} ADU")
            else:
                print(f"  ✓ Results match within tolerance")

        # Save diagnostic output
        diagnostic_file = tmp_path / "diagnostic.txt"
        with open(diagnostic_file, "w") as f:
            f.write("="*80 + "\n")
            f.write("BROADENED PSF DIAGNOSTIC\n")
            f.write("="*80 + "\n")
            f.write(f"Template shape: {template.shape}\n")
            f.write(f"Science shape: {science.shape}\n")
            f.write(f"Template PSF σ: {PSF_SIGMA_TEMPLATE} px\n")
            f.write(f"Science PSF σ: {PSF_SIGMA_SCIENCE:.3f} px\n")
            f.write(f"Stars: {len(star_field['fluxes'])}\n")
            f.write(f"\nKERNEL FIT RESULTS:\n")
            f.write(f"  Kernel norm: {kernel_solution.kernel_norm:.6f}\n")
            f.write(f"  Valid stamps (nS): {nS}\n")
            f.write(f"  Mean sigma: {kernel_solution.mean_sigma:.6f}\n")
            f.write(f"  Scatter sigma: {kernel_solution.scatter_sigma:.6f}\n")
            f.write(f"\nDIFFERENCE IMAGE STATS:\n")
            f.write(f"  CLI  : min={cli_diff.min():.2f}, max={cli_diff.max():.2f}, "
                   f"mean={cli_diff.mean():.2f}\n")
            f.write(f"  API  : min={api_diff.min():.2f}, max={api_diff.max():.2f}, "
                   f"mean={api_diff.mean():.2f}\n")
            f.write(f"  Diff : min={pixel_diffs.min():.6f}, max={pixel_diffs.max():.6f}, "
                   f"mean={pixel_diffs.mean():.6f}\n")

        print(f"\nDiagnostic output saved to: {diagnostic_file}")

        # Assertion: should match within tolerance
        np.testing.assert_allclose(
            cli_diff, api_diff,
            rtol=1e-3, atol=0.5,
            err_msg=f"Python API and CLI produce different results (max diff {pixel_diffs.max():.2f} ADU)"
        )
