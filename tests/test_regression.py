"""
End-to-end regression tests for HOTPANTS.

Each test builds a well-understood synthetic scenario, runs hotpants, and
checks that the output satisfies a statistical invariant.  The goal is to
document *current* behaviour so that regressions are caught when the C code
is refactored or its dependencies are replaced.

Coordinate convention: numpy arrays are (row=y, col=x).
"""

import numpy as np
import pytest

from .conftest import (
    BACKGROUND,
    PSF_SIGMA_EXTRA,
    PSF_SIGMA_SCIENCE,
    PSF_SIGMA_TEMPLATE,
    RDNOISE,
    NX,
    NY,
    load_diff,
    make_image,
    run_hotpants,
    write_fits,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def expected_noise(n_images=2):
    """Rough per-pixel noise floor in a difference of n_images synthetic frames."""
    variance_per_pixel = BACKGROUND + RDNOISE**2  # dominant at faint sky level
    return np.sqrt(n_images * variance_per_pixel)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestIdenticalImages:
    """
    Template and science are different noise realisations of the same scene.

    After subtracting, the residual should be pure noise — no coherent
    signal, variance consistent with 2× the per-pixel noise of the inputs.
    """

    def test_returncode(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._run(hotpants_binary, star_field, tmp_path)
        result = run_hotpants(hotpants_binary, tmpl, sci, diff)
        assert result.returncode == 0, result.stderr

    def test_mean_consistent_with_zero(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._run(hotpants_binary, star_field, tmp_path)
        run_hotpants(hotpants_binary, tmpl, sci, diff)
        d = load_diff(diff)
        # Mean should not exceed 5-sigma from zero
        assert abs(np.mean(d)) < 5 * np.std(d) / np.sqrt(len(d))

    def test_noise_level(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._run(hotpants_binary, star_field, tmp_path)
        run_hotpants(hotpants_binary, tmpl, sci, diff)
        d = load_diff(diff)
        # Std should be within a factor of 3 of expected noise
        assert np.std(d) < 3 * expected_noise(2)

    @staticmethod
    def _run(binary, star_field, tmp_path):
        tmpl = tmp_path / "tmpl.fits"
        sci = tmp_path / "sci.fits"
        diff = tmp_path / "diff.fits"
        noiseless = star_field["noiseless"]
        write_fits(tmpl, make_image(noiseless, PSF_SIGMA_TEMPLATE, np.random.default_rng(1)))
        write_fits(sci,  make_image(noiseless, PSF_SIGMA_TEMPLATE, np.random.default_rng(2)))
        return tmpl, sci, diff


class TestGaussianKernelRecovery:
    """
    Science PSF is template PSF broadened by an additional Gaussian.

    Hotpants must discover and apply the connecting kernel so that the
    difference image collapses to noise.
    """

    def test_returncode(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        result = run_hotpants(hotpants_binary, tmpl, sci, diff)
        assert result.returncode == 0, result.stderr

    def test_difference_consistent_with_noise(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        run_hotpants(hotpants_binary, tmpl, sci, diff)
        d = load_diff(diff)
        # Residual std should be within a factor of 5 of the noise floor;
        # factor is generous to tolerate fitting artefacts near bright stars.
        assert np.std(d) < 5 * expected_noise(2)

    def test_no_large_residuals(self, hotpants_binary, star_field, tmp_path):
        """No pixel should deviate by more than 20-sigma from the mean residual."""
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        run_hotpants(hotpants_binary, tmpl, sci, diff)
        d = load_diff(diff)
        sigma = np.std(d)
        assert np.max(np.abs(d - np.mean(d))) < 20 * sigma

    @staticmethod
    def _inputs(star_field, tmp_path):
        from scipy.ndimage import gaussian_filter
        tmpl = tmp_path / "tmpl.fits"
        sci = tmp_path / "sci.fits"
        diff = tmp_path / "diff.fits"
        noiseless = star_field["noiseless"]
        # Template: template PSF only
        write_fits(tmpl, make_image(noiseless, PSF_SIGMA_TEMPLATE, np.random.default_rng(10)))
        # Science: convolve the noiseless field with the extra broadening kernel
        # first, then add the template PSF on top, giving effective sigma = PSF_SIGMA_SCIENCE
        from scipy.ndimage import gaussian_filter as gf
        pre_broadened = gf(noiseless, sigma=PSF_SIGMA_EXTRA, mode="constant")
        write_fits(sci, make_image(pre_broadened, PSF_SIGMA_TEMPLATE, np.random.default_rng(11)))
        return tmpl, sci, diff


class TestPointSourceInjection:
    """
    A synthetic transient is present in the science image but not the template.

    After subtraction the transient should appear as a significant peak in
    the difference image at its known position.
    """

    TRANSIENT_X = NX // 2
    TRANSIENT_Y = NY // 2
    TRANSIENT_FLUX = 25_000.0
    SEARCH_RADIUS = 12  # pixels around injected position

    def test_returncode(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        result = run_hotpants(hotpants_binary, tmpl, sci, diff)
        assert result.returncode == 0, result.stderr

    def test_transient_detected(self, hotpants_binary, star_field, tmp_path):
        """Peak near the transient position must be > 5σ above the background noise."""
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        run_hotpants(hotpants_binary, tmpl, sci, diff)
        d = load_diff(diff)
        noise = np.std(d)

        from astropy.io import fits
        diff_map = fits.getdata(str(diff)).astype(np.float64)
        r = self.SEARCH_RADIUS
        ty, tx = self.TRANSIENT_Y, self.TRANSIENT_X
        region = diff_map[ty - r : ty + r, tx - r : tx + r]
        assert np.nanmax(region) > 5 * noise

    def test_transient_is_local_maximum(self, hotpants_binary, star_field, tmp_path):
        """The transient region peak should be the single highest pixel in the diff."""
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        run_hotpants(hotpants_binary, tmpl, sci, diff)

        from astropy.io import fits
        diff_map = fits.getdata(str(diff)).astype(np.float64)
        r = self.SEARCH_RADIUS
        ty, tx = self.TRANSIENT_Y, self.TRANSIENT_X
        region = diff_map[ty - r : ty + r, tx - r : tx + r]
        outer = diff_map.copy()
        outer[ty - r : ty + r, tx - r : tx + r] = np.nan
        assert np.nanmax(region) > np.nanmax(outer)

    @classmethod
    def _inputs(cls, star_field, tmp_path):
        tmpl = tmp_path / "tmpl.fits"
        sci = tmp_path / "sci.fits"
        diff = tmp_path / "diff.fits"
        noiseless = star_field["noiseless"]
        write_fits(tmpl, make_image(noiseless, PSF_SIGMA_TEMPLATE, np.random.default_rng(20)))
        write_fits(
            sci,
            make_image(
                noiseless,
                PSF_SIGMA_TEMPLATE,
                np.random.default_rng(21),
                extra_sources=[(cls.TRANSIENT_X, cls.TRANSIENT_Y, cls.TRANSIENT_FLUX)],
            ),
        )
        return tmpl, sci, diff


class TestPhotometricNormalization:
    """
    Science image has a different (known) flux scale from the template.

    After hotpants normalises to the template photometric system the
    difference should collapse to noise.
    """

    SCALE = 1.5  # science flux = 1.5 × template flux

    def test_returncode(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        result = run_hotpants(
            hotpants_binary, tmpl, sci, diff, extra_args=["-n", "t"]
        )
        assert result.returncode == 0, result.stderr

    def test_difference_consistent_with_noise(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        run_hotpants(hotpants_binary, tmpl, sci, diff, extra_args=["-n", "t"])
        d = load_diff(diff)
        assert np.std(d) < 5 * expected_noise(2)

    @classmethod
    def _inputs(cls, star_field, tmp_path):
        tmpl = tmp_path / "tmpl.fits"
        sci = tmp_path / "sci.fits"
        diff = tmp_path / "diff.fits"
        noiseless = star_field["noiseless"]
        write_fits(tmpl, make_image(noiseless, PSF_SIGMA_TEMPLATE, np.random.default_rng(30)))
        write_fits(
            sci,
            make_image(
                noiseless,
                PSF_SIGMA_TEMPLATE,
                np.random.default_rng(31),
                scale=cls.SCALE,
            ),
        )
        return tmpl, sci, diff


class TestMultiRegion:
    """
    Hotpants can split the image into spatial regions and fit an independent
    kernel per region.  With a 2×2 split each region contains ~20 stars,
    which is enough to fit a spatially constant kernel (order 0) but not a
    higher-order polynomial.  The test deliberately uses a simple same-PSF
    scenario to isolate region-splitting from PSF-matching.
    """

    # Reduce spatial kernel order to 0 (constant) so the fit is fully
    # determined even with the small number of stamps per region.
    REGION_ARGS = ["-nrx", "2", "-nry", "2", "-ko", "0", "-bgo", "0"]

    def test_returncode(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        result = run_hotpants(hotpants_binary, tmpl, sci, diff, extra_args=self.REGION_ARGS)
        assert result.returncode == 0, result.stderr

    def test_difference_consistent_with_noise(self, hotpants_binary, star_field, tmp_path):
        tmpl, sci, diff = self._inputs(star_field, tmp_path)
        run_hotpants(hotpants_binary, tmpl, sci, diff, extra_args=self.REGION_ARGS)
        d = load_diff(diff)
        assert np.std(d) < 5 * expected_noise(2)

    @staticmethod
    def _inputs(star_field, tmp_path):
        tmpl = tmp_path / "tmpl.fits"
        sci = tmp_path / "sci.fits"
        diff = tmp_path / "diff.fits"
        noiseless = star_field["noiseless"]
        # Use same PSF for both images — the test is about region splitting,
        # not kernel estimation accuracy.
        write_fits(tmpl, make_image(noiseless, PSF_SIGMA_TEMPLATE, np.random.default_rng(40)))
        write_fits(sci,  make_image(noiseless, PSF_SIGMA_TEMPLATE, np.random.default_rng(41)))
        return tmpl, sci, diff
