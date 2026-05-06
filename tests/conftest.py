import subprocess
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from scipy.ndimage import gaussian_filter

REPO_ROOT = Path(__file__).parent.parent
HOTPANTS_BIN = REPO_ROOT / "build" / "hotpants"

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
    """Return path to the hotpants binary, building it via CMake if necessary."""
    if not HOTPANTS_BIN.exists():
        result = subprocess.run(
            ["cmake", "-B", "build", "-DCMAKE_BUILD_TYPE=Release"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            pytest.fail(f"cmake configure failed:\n{result.stderr}")
        result = subprocess.run(
            ["cmake", "--build", "build", "--parallel"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            pytest.fail(f"cmake build failed:\n{result.stderr}")
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
    """Run hotpants; raise error if command fails."""
    cmd = (
        [str(binary)]
        + ["-inim", str(sci), "-tmplim", str(tmpl), "-outim", str(diff)]
        + HOTPANTS_BASE_ARGS
        + list(extra_args)
    )
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        msg = f"hotpants failed with exit code {result.returncode}\nStdout: {result.stdout}\nStderr: {result.stderr}"
        raise RuntimeError(msg)
    return result


def load_diff(path):
    """Load difference image, masking hotpants fill pixels (1e-30)."""
    data = fits.getdata(str(path)).astype(np.float64)
    return data[np.abs(data - 1e-30) > 1e-35]


# ---------------------------------------------------------------------------
# Python API Test Fixtures (Mock C Library)
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def mock_hotpants_library(monkeypatch, request):
    """
    Auto-mock the HOTPANTS C library for Python API unit tests.

    Provides dummy implementations of C library functions so tests can run
    without requiring the compiled C library to be available.

    This fixture is automatically applied to all tests in test_api.py.
    Integration tests (test_api_integration.py) use the real C library.
    """
    # Only mock for test_api.py (unit tests), not integration or regression tests
    if 'test_api_integration' in str(request.fspath) or 'test_regression' in str(request.fspath):
        return None
    if 'test_api.py' not in str(request.fspath):
        return None

    from unittest.mock import MagicMock
    from hotpants import _hotpants_ffi

    # Mock the C library loader to return a mock library
    mock_lib = MagicMock()

    # Mock global variable getters/setters
    # These are called by _core.py when setting up the global state
    globals_state = {
        'hwKernel': 15,
        'kerOrder': 2,
        'bgOrder': 1,
        'nKSStamps': 3,
        'hwKSStamp': 10,
        'nRegX': 1,
        'nRegY': 1,
        'nStampX': 10,
        'nStampY': 10,
        'useFullSS': 0,
        'tUThresh': 25000.0,
        'tLThresh': 0.0,
        'iUThresh': 25000.0,
        'iLThresh': 0.0,
        'tGain': 1.0,
        'iGain': 1.0,
        'tRdnoise': 0.0,
        'iRdnoise': 0.0,
        'tPedestal': 0.0,
        'iPedestal': 0.0,
        'nCompKer': 3,
        'nComp': 6,  # (2+1)*(2+2)/2 for kerOrder=2
        'nCompBG': 3,  # (1+1)*(1+2)/2 for bgOrder=1
        'verbose': 0,
        'nThread': 1,
    }

    def mock_get_global_int(name):
        return globals_state.get(name, 0)

    def mock_set_global_int(name, value):
        globals_state[name] = value

    def mock_get_global_float(name):
        return float(globals_state.get(name, 0.0))

    def mock_set_global_float(name, value):
        globals_state[name] = value

    # Mock C function return values for kernel fitting
    # allocateStamps should return 0 (success)
    mock_lib.allocateStamps.return_value = 0
    # fitKernel returns void (None)
    mock_lib.fitKernel.return_value = None
    # spatial_convolve returns void (None)
    mock_lib.spatial_convolve.return_value = None
    # buildStamps returns void (None)
    mock_lib.buildStamps.return_value = None
    # freeStampMem returns void (None)
    mock_lib.freeStampMem.return_value = None
    # Wrapper functions
    # initBuildStampsContext should return 0 (success)
    mock_lib.initBuildStampsContext.return_value = 0
    # buildStampsRegion should return 0 (success) and set output pointers
    mock_lib.buildStampsRegion.return_value = 0
    # buildStamps returns void (None)
    mock_lib.buildStamps.return_value = None
    # cleanupBuildStampsContext returns void (None)
    mock_lib.cleanupBuildStampsContext.return_value = None

    # Patch the library loading and global variable functions
    monkeypatch.setattr(_hotpants_ffi, 'get_library', lambda: mock_lib)
    monkeypatch.setattr(_hotpants_ffi, 'get_global_int', mock_get_global_int)
    monkeypatch.setattr(_hotpants_ffi, 'set_global_int', mock_set_global_int)
    monkeypatch.setattr(_hotpants_ffi, 'get_global_float', mock_get_global_float)
    monkeypatch.setattr(_hotpants_ffi, 'set_global_float', mock_set_global_float)

    return mock_lib
