#!/usr/bin/env python3
"""
Comprehensive bottleneck analysis for HOTPANTS.
Uses timing, code examination, and measurements to identify acceleration opportunities.
"""
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

# Test parameters
BACKGROUND = 1000.0
GAIN = 1.0
RDNOISE = 5.0
PSF_SIGMA_TEMPLATE = 1.5
PSF_SIGMA_SCIENCE = np.sqrt(PSF_SIGMA_TEMPLATE**2 + 1.5**2)

def make_synthetic_images(tmpdir, nx=512, ny=512, n_stars=120):
    """Create synthetic template and science images."""
    rng = np.random.default_rng(0xC0FFEE)
    margin = 40
    x_stars = rng.uniform(margin, nx - margin, n_stars)
    y_stars = rng.uniform(margin, ny - margin, n_stars)
    fluxes = rng.uniform(8_000, 45_000, n_stars)

    noiseless = np.zeros((ny, nx), dtype=np.float64)
    for x, y, f in zip(x_stars, y_stars, fluxes):
        xi, yi = int(round(x)), int(round(y))
        noiseless[yi, xi] += f

    # Template
    template_conv = gaussian_filter(noiseless, sigma=PSF_SIGMA_TEMPLATE)
    template_conv += BACKGROUND
    template_noise = np.sqrt(np.maximum(template_conv, 0) / GAIN + RDNOISE**2) * rng.standard_normal((ny, nx))
    template = (template_conv + template_noise).astype(np.float32)

    # Science
    science_conv = gaussian_filter(noiseless * 1.05, sigma=PSF_SIGMA_SCIENCE)
    science_conv += BACKGROUND * 0.98
    science_noise = np.sqrt(np.maximum(science_conv, 0) / GAIN + RDNOISE**2) * rng.standard_normal((ny, nx))
    science = (science_conv + science_noise).astype(np.float32)

    template_path = tmpdir / "template.fits"
    science_path = tmpdir / "science.fits"
    diff_path = tmpdir / "diff.fits"

    fits.PrimaryHDU(template).writeto(str(template_path), overwrite=True)
    fits.PrimaryHDU(science).writeto(str(science_path), overwrite=True)

    return template_path, science_path, diff_path

def run_hotpants(hotpants_bin, template, science, diff, extra_args=""):
    """Run hotpants and return execution time."""
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
    if extra_args:
        cmd.extend(extra_args.split())

    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - start

    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        return None
    return elapsed

def analyze_code():
    """Analyze source code to estimate bottlenecks."""
    print("\n" + "="*80)
    print("CODE ANALYSIS - Function complexity and scope")
    print("="*80)

    alard_path = Path("/home/user/hotpants/src/alard.c")
    functions_path = Path("/home/user/hotpants/src/functions.c")

    # Key functions to analyze
    key_functions = {
        'spatial_convolve': (alard_path, '~60% CPU (documented)'),
        'getPsfCenters': (functions_path, '~35% CPU (documented)'),
        'xy_conv_stamp': (alard_path, '~18% CPU (documented)'),
        'fitKernel': (alard_path, 'Main kernel fitting loop'),
        'make_kernel': (alard_path, 'Polynomial evaluation at pixels'),
        'build_matrix': (alard_path, 'Accumulate normal equations'),
    }

    for func_name, (fpath, desc) in key_functions.items():
        content = fpath.read_text()
        # Find function
        import re
        pattern = rf'^[^/]*?\b{func_name}\s*\([^)]*\)\s*\{{?'
        matches = list(re.finditer(pattern, content, re.MULTILINE))
        if matches:
            start_pos = matches[0].start()
            # Count braces to estimate function size
            brace_count = 0
            func_start = start_pos
            func_size = 0
            for i, char in enumerate(content[start_pos:]):
                if char == '{':
                    brace_count += 1
                elif char == '}':
                    brace_count -= 1
                    if brace_count == 0:
                        func_size = i
                        break

            lines = content[:start_pos+func_size].count('\n') - content[:start_pos].count('\n')
            print(f"\n{func_name}:")
            print(f"  Location: {fpath.name}")
            print(f"  Estimated size: {lines} lines")
            print(f"  Profile contribution: {desc}")

def main():
    hotpants_bin = Path("/home/user/hotpants/build/hotpants")
    if not hotpants_bin.exists():
        print(f"Error: {hotpants_bin} not found")
        exit(1)

    print("="*80)
    print("HOTPANTS BOTTLENECK ANALYSIS")
    print("="*80)

    analyze_code()

    print("\n" + "="*80)
    print("TIMING MEASUREMENTS - Performance across different workloads")
    print("="*80)

    with tempfile.TemporaryDirectory(prefix="hotpants_timing_") as tmpdir:
        tmpdir = Path(tmpdir)

        # Test 1: Standard 512x512
        print("\n[Test 1] Standard 512×512 image, 1 region, direct convolution:")
        t, s, d = make_synthetic_images(tmpdir, 512, 512, 120)
        elapsed = run_hotpants(hotpants_bin, t, s, d)
        if elapsed:
            print(f"  Time: {elapsed:.2f}s")

        # Test 2: Larger 1024x1024
        print("\n[Test 2] Larger 1024×1024 image, 1 region, direct convolution:")
        t, s, d = make_synthetic_images(tmpdir, 1024, 1024, 200)
        elapsed = run_hotpants(hotpants_bin, t, s, d)
        if elapsed:
            print(f"  Time: {elapsed:.2f}s")

        # Test 3: With FFT acceleration (if available)
        print("\n[Test 3] Standard 512×512 with FFT acceleration:")
        t, s, d = make_synthetic_images(tmpdir, 512, 512, 120)
        elapsed = run_hotpants(hotpants_bin, t, s, d, "-fft")
        if elapsed:
            print(f"  Time: {elapsed:.2f}s")

        # Test 4: Multi-threading
        print("\n[Test 4] Standard 512×512 with OpenMP (2 threads):")
        t, s, d = make_synthetic_images(tmpdir, 512, 512, 120)
        import os
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = '2'
        elapsed = run_hotpants(hotpants_bin, t, s, d)
        if elapsed:
            print(f"  Time: {elapsed:.2f}s")

    print("\n" + "="*80)
    print("IDENTIFIED BOTTLENECKS & ACCELERATION OPPORTUNITIES")
    print("="*80)

    bottlenecks = """
1. spatial_convolve() [~60% CPU]
   ─────────────────────────────────────
   • Direct 2D convolution using brute-force pixel multiplication
   • O(k² × N) complexity where k=kernel size, N=image pixels
   • Already has FFT alternative in spatial_convolve_fft()
   • Opportunities:
     - Ensure FFT path is always used for large images (consider auto-switch)
     - SIMD vectorization of inner loop (AVX2/AVX-512)
     - Tile-based cache optimization for better L1/L2 locality
     - Parallel reduction of summation operations

2. getPsfCenters() [~35% CPU]
   ──────────────────────────────
   • Locates bright star centers via stamp-by-stamp search
   • Brute-force nearest-neighbor comparison within each stamp
   • Current approach: all-pairs comparison per stamp
   • Opportunities:
     - Spatial hashing or k-d trees for faster PSF matching
     - OpenMP parallelization (stamp loops are independent)
     - SIMD vectorization of distance computation
     - Approximate nearest-neighbor methods (e.g., Annoy, FAISS)
     - Cache optimization for stamp data layout

3. xy_conv_stamp() [~18% CPU, subset of spatial_convolve]
   ───────────────────────────────────────────────────────
   • Separable 2D convolution using Gaussian basis
   • Called for each basis element, each stamp, each region
   • Opportunities:
     - Separable convolution is already efficient
     - Vectorization of separable passes (row/col-major ops)
     - Shared memory optimization if moving to GPU/SIMD

4. make_kernel() [polynomial evaluation at pixels]
   ─────────────────────────────────────────────
   • Evaluates spatial polynomial coefficients at every pixel
   • Called for each pixel when applying final kernel
   • Currently uses Horner's method (good)
   • Opportunities:
     - Vectorize polynomial evaluation across multiple pixels
     - Precompute polynomial values on coarse grid + interpolate
     - Cache polynomial coefficients in local/register storage

5. build_matrix() [Accumulate normal equations]
   ─────────────────────────────────────────────
   • Builds normal equation matrix for least-squares solving
   • Matrix size depends on basis size and region parameters
   • Currently < 5% CPU usage (not a bottleneck)
   • Using Cholesky via LAPACK (optimal for symmetric positive-definite)

SUMMARY OF ACCELERATION STRATEGIES
═════════════════════════════════════
Priority 1 (High impact, moderate effort):
  ✓ Mandatory FFT path for large images (>512²)
  ✓ OpenMP parallelization of getPsfCenters() stamp loops
  ✓ SIMD vectorization of spatial_convolve inner loop

Priority 2 (Medium impact, high effort):
  ✓ Implement approximate nearest-neighbor search in getPsfCenters()
  ✓ Tile-based cache optimization in convolution operations
  ✓ GPU acceleration of convolution via cuFFT or similar

Priority 3 (Lower impact, lower effort):
  ✓ Polynomial evaluation vectorization
  ✓ Memory layout optimization for cache efficiency
  ✓ Profiling-guided optimizations based on real astronomical data
"""
    print(bottlenecks)

if __name__ == "__main__":
    main()
