#!/usr/bin/env python3
"""
Profiling script for HOTPANTS using gprof.
Generates synthetic FITS images and runs hotpants with gprof profiling.
"""
import subprocess
import tempfile
from pathlib import Path
import os

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

def profile_hotpants(hotpants_bin, template, science, diff, workdir):
    """Run hotpants with gprof for profiling."""
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

    print(f"Running: {' '.join(cmd)}")
    env = os.environ.copy()
    env['GMON_OUT_PREFIX'] = str(workdir / 'gmon')
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(workdir), env=env)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        return False

    print(f"Hotpants execution completed")
    return True

def analyze_gprof(workdir):
    """Analyze the gprof profile."""
    gmon_file = workdir / 'gmon.out'
    if not gmon_file.exists():
        print(f"Error: gmon.out not found in {workdir}")
        return

    hotpants_bin = Path("/home/user/hotpants/build/hotpants")

    # Run gprof
    result = subprocess.run(
        ["gprof", str(hotpants_bin), str(gmon_file)],
        capture_output=True,
        text=True
    )

    print("\n" + "="*80)
    print("GPROF PROFILE (Flat profile - Top 40 functions by CPU time)")
    print("="*80)

    # Parse and format output
    lines = result.stdout.split('\n')
    in_flat = False
    count = 0
    for i, line in enumerate(lines):
        if 'Flat profile' in line:
            in_flat = True
            print(line)
            continue

        if in_flat:
            if line.startswith('Call graph') or line.startswith('['):
                break
            if 'cumulative' in line or '%' in line or line.strip() == '':
                print(line)
                count += 1
                if count > 45:
                    break

    print("\n" + "="*80)
    print("GPROF CALL GRAPH (Top 60 lines)")
    print("="*80)

    in_graph = False
    count = 0
    for line in lines:
        if 'Call graph' in line:
            in_graph = True
        if in_graph:
            print(line)
            count += 1
            if count > 60:
                break

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
        print(f"  Diff:     {diff}\n")

        if profile_hotpants(hotpants_bin, template, science, diff, tmpdir):
            analyze_gprof(tmpdir)

        print("\nDone!")
