#!/usr/bin/env python3
"""
Profiling script for HOTPANTS to identify bottlenecks.
Generates synthetic FITS images and runs hotpants with profiling.
"""
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

# Test image parameters
NY, NX = 512, 512
BACKGROUND = 1000.0
GAIN = 1.0
RDNOISE = 5.0
PSF_SIGMA_TEMPLATE = 1.5
PSF_SIGMA_SCIENCE = np.sqrt(PSF_SIGMA_TEMPLATE**2 + 1.5**2)

def make_synthetic_images(tmpdir):
    """Create synthetic template and science images for profiling."""
    rng = np.random.default_rng(0xC0FFEE)
    n_stars = 120
    margin = 40
    x_stars = rng.uniform(margin, NX - margin, n_stars)
    y_stars = rng.uniform(margin, NY - margin, n_stars)
    fluxes = rng.uniform(8_000, 45_000, n_stars)

    # Create point source grid
    noiseless = np.zeros((NY, NX), dtype=np.float64)
    for x, y, f in zip(x_stars, y_stars, fluxes):
        xi, yi = int(round(x)), int(round(y))
        noiseless[yi, xi] += f

    # Template image (narrower PSF)
    template_conv = gaussian_filter(noiseless, sigma=PSF_SIGMA_TEMPLATE)
    template_conv += BACKGROUND
    template_noise = np.sqrt(np.maximum(template_conv, 0) / GAIN + RDNOISE**2) * rng.standard_normal((NY, NX))
    template = (template_conv + template_noise).astype(np.float32)

    # Science image (broader PSF, scaled slightly differently)
    science_conv = gaussian_filter(noiseless * 1.05, sigma=PSF_SIGMA_SCIENCE)
    science_conv += BACKGROUND * 0.98
    science_noise = np.sqrt(np.maximum(science_conv, 0) / GAIN + RDNOISE**2) * rng.standard_normal((NY, NX))
    science = (science_conv + science_noise).astype(np.float32)

    # Write FITS files
    template_path = tmpdir / "template.fits"
    science_path = tmpdir / "science.fits"
    diff_path = tmpdir / "diff.fits"

    fits.PrimaryHDU(template).writeto(str(template_path), overwrite=True)
    fits.PrimaryHDU(science).writeto(str(science_path), overwrite=True)

    return template_path, science_path, diff_path

def profile_hotpants(hotpants_bin, template, science, diff):
    """Run hotpants with perf for profiling."""
    cmd = [
        str(hotpants_bin),
        "-inim", str(science),
        "-tmplim", str(template),
        "-outim", str(diff),
        "-tu", "60000", "-tl", "-200",
        "-iu", "60000", "-il", "-200",
        "-tg", str(GAIN), "-ig", str(GAIN),
        "-tr", str(RDNOISE), "-ir", str(RDNOISE),
        "-r", "8",
        "-nsx", "5", "-nsy", "5",
        "-nss", "3",
        "-ft", "8",
        "-v", "0",
    ]

    # Run with perf record
    perf_cmd = ["perf", "record", "-g", "-o", "/tmp/perf.data"] + cmd
    print(f"Running: {' '.join(perf_cmd)}")
    result = subprocess.run(perf_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        return False

    print("\nHotpants execution completed. Analyzing profile...")
    return True

def analyze_profile():
    """Analyze the perf profile."""
    # Generate report
    result = subprocess.run(
        ["perf", "report", "-i", "/tmp/perf.data", "--stdio", "--no-children"],
        capture_output=True,
        text=True
    )

    print("\n" + "="*80)
    print("PERF REPORT (Top functions by CPU time)")
    print("="*80)
    print(result.stdout[:3000])  # Print first 3000 chars of output

    # Generate flame graph data if possible
    print("\n" + "="*80)
    print("Attempting to generate call graph...")
    print("="*80)
    result = subprocess.run(
        ["perf", "report", "-i", "/tmp/perf.data", "-g", "caller", "--stdio"],
        capture_output=True,
        text=True,
        timeout=10
    )
    print(result.stdout[:2000])

if __name__ == "__main__":
    hotpants_bin = Path("/home/user/hotpants/build/hotpants")
    if not hotpants_bin.exists():
        print(f"Error: {hotpants_bin} not found. Build the project first.")
        exit(1)

    with tempfile.TemporaryDirectory(prefix="hotpants_profile_") as tmpdir:
        tmpdir = Path(tmpdir)
        print("Generating synthetic FITS images...")
        template, science, diff = make_synthetic_images(tmpdir)
        print(f"  Template: {template}")
        print(f"  Science:  {science}")
        print(f"  Diff:     {diff}")

        if profile_hotpants(hotpants_bin, template, science, diff):
            analyze_profile()

        print("\nDone! Profile data saved to /tmp/perf.data")
