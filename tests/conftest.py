import subprocess
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from scipy.ndimage import gaussian_filter

REPO_ROOT = Path(__file__).parent.parent
HOTPANTS_BIN = REPO_ROOT / "hotpants"

# Fixed image geometry for all synthetic tests
NY, NX = 512, 512
BACKGROUND = 1000.0
GAIN = 1.0
RDNOISE = 5.0

# PSF widths: science image has broader PSF than template
PSF_SIGMA_TEMPLATE = 1.5   # pixels
PSF_SIGMA_EXTRA = 1.5      # kernel connecting template to science
PSF_SIGMA_SCIENCE = np.sqrt(PSF_SIGMA_TEMPLATE**2 + PSF_SIGMA_EXTRA**2)

# Default hotpants options for test images
HOTPANTS_BASE_ARGS = [
    # Use thresholds suited to synthetic images
    "-tu", "60000", "-tl", "-200",
    "-iu", "60000", "-il", "-200",
    # Calibration
    "-tg", str(GAIN), "-ig", str(GAIN),
    "-tr", str(RDNOISE), "-ir", str(RDNOISE),
    # Kernel geometry (smaller than default for speed)
    "-r", "8",
    # Stamp grid: 5×5 per region for 512×512 image
    "-nsx", "5", "-nsy", "5",
    # Number of substamps per stamp
    "-nss", "3",
    # Lower fit threshold to work with synthetic data
    "-ft", "8",
    # Suppress most output
    "-v", "0",
]


@pytest.fixture(scope="session")
def hotpants_binary():
    """Return path to the hotpants binary, building it first if necessary."""
    result = subprocess.run(
        ["make", "hotpants"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(f"hotpants build failed:\n{result.stderr}")
    if not HOTPANTS_BIN.exists():
        pytest.fail("hotpants binary not found after build")
    return HOTPANTS_BIN


@pytest.fixture(scope="session")
def star_field():
    """
    Reproducible synthetic star field.

    Returns a dict with the noiseless flux image (point sources, pre-convolution)
    and the star catalogue.
    """
    rng = np.random.default_rng(0xC0FFEE)
    n_stars = 80
    # Keep stars away from borders so PSF convolution has room
    margin = 40
    x_stars = rng.uniform(margin, NX - margin, n_stars)
    y_stars = rng.uniform(margin, NY - margin, n_stars)
    # Flux range chosen so even the faintest star is well above the fit threshold
    fluxes = rng.uniform(8_000, 45_000, n_stars)

    # Place delta-function sources on the grid
    noiseless = np.zeros((NY, NX), dtype=np.float64)
    for x, y, f in zip(x_stars, y_stars, fluxes):
        xi, yi = int(round(x)), int(round(y))
        noiseless[yi, xi] += f

    return {
        "noiseless": noiseless,
        "x_stars": x_stars,
        "y_stars": y_stars,
        "fluxes": fluxes,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_image(noiseless, psf_sigma, rng, *, scale=1.0, extra_sources=()):
    """
    Build a realistic simulated FITS image.

    Parameters
    ----------
    noiseless : ndarray
        Point-source flux map (delta functions, pre-PSF).
    psf_sigma : float
        Gaussian PSF sigma in pixels.
    rng : numpy.random.Generator
    scale : float
        Multiplicative flux scale applied before noise (to simulate
        photometric differences).
    extra_sources : iterable of (x, y, flux)
        Additional point sources added before convolution (e.g. transients).
    """
    src = noiseless.copy()
    for x, y, flux in extra_sources:
        xi, yi = int(round(x)), int(round(y))
        src[yi, xi] += flux

    # Convolve with PSF
    convolved = gaussian_filter(src, sigma=psf_sigma, mode="constant") * scale
    # Add sky background
    convolved += BACKGROUND
    # Poisson noise from signal + sky, plus Gaussian read noise
    variance = np.maximum(convolved, 0) / GAIN + RDNOISE**2
    noise = np.sqrt(variance) * rng.standard_normal(convolved.shape)
    return (convolved + noise).astype(np.float32)


def write_fits(path, data):
    fits.PrimaryHDU(data).writeto(str(path), overwrite=True)


def run_hotpants(binary, tmpl, sci, diff, extra_args=()):
    """Run hotpants; return the CompletedProcess."""
    cmd = (
        [str(binary)]
        + ["-inim", str(sci), "-tmplim", str(tmpl), "-outim", str(diff)]
        + HOTPANTS_BASE_ARGS
        + list(extra_args)
    )
    return subprocess.run(cmd, capture_output=True, text=True)


def load_diff(path):
    """Load difference image, masking hotpants fill pixels (1e-30)."""
    data = fits.getdata(str(path)).astype(np.float64)
    return data[np.abs(data - 1e-30) > 1e-35]
