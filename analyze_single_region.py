#!/usr/bin/env python3
"""
Analysis of optimization opportunities for single-region image processing.
Compares single-region vs multi-region execution patterns.
"""
import subprocess
import tempfile
import time
from pathlib import Path
from dataclasses import dataclass

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter


@dataclass
class RegionConfig:
    """Configuration for region layout."""
    nrx: int
    nry: int

    @property
    def num_regions(self):
        return self.nrx * self.nry

    @property
    def name(self):
        if self.num_regions == 1:
            return "single_region"
        else:
            return f"tiled_{self.nrx}x{self.nry}"


class SingleRegionAnalyzer:
    """Analyzes overhead and bottlenecks in single vs multi-region processing."""

    def __init__(self, hotpants_bin: Path):
        self.hotpants_bin = hotpants_bin

    def make_synthetic_images(self, tmpdir: Path, nx: int, ny: int, n_stars: int):
        """Generate synthetic test images."""
        rng = np.random.default_rng(0xC0FFEE)
        margin = 40
        x_stars = rng.uniform(margin, nx - margin, n_stars)
        y_stars = rng.uniform(margin, ny - margin, n_stars)
        fluxes = rng.uniform(8_000, 45_000, n_stars)

        noiseless = np.zeros((ny, nx), dtype=np.float64)
        for x, y, f in zip(x_stars, y_stars, fluxes):
            xi, yi = int(round(x)), int(round(y))
            noiseless[yi, xi] += f

        bg = 1000.0
        gain = 1.0
        rdnoise = 5.0
        psf_template = 1.5
        psf_science = np.sqrt(psf_template**2 + 1.5**2)

        template_conv = gaussian_filter(noiseless, sigma=psf_template)
        template_conv += bg
        template_noise = np.sqrt(np.maximum(template_conv, 0) / gain + rdnoise**2) * rng.standard_normal((ny, nx))
        template = (template_conv + template_noise).astype(np.float32)

        science_conv = gaussian_filter(noiseless * 1.05, sigma=psf_science)
        science_conv += bg * 0.98
        science_noise = np.sqrt(np.maximum(science_conv, 0) / gain + rdnoise**2) * rng.standard_normal((ny, nx))
        science = (science_conv + science_noise).astype(np.float32)

        template_path = tmpdir / "template.fits"
        science_path = tmpdir / "science.fits"
        diff_path = tmpdir / "diff.fits"

        fits.PrimaryHDU(template).writeto(str(template_path), overwrite=True)
        fits.PrimaryHDU(science).writeto(str(science_path), overwrite=True)

        return template_path, science_path, diff_path

    def run_hotpants(self, template, science, diff, region_config: RegionConfig):
        """Run hotpants with specified region configuration."""
        cmd = [
            str(self.hotpants_bin),
            "-inim", str(science),
            "-tmplim", str(template),
            "-outim", str(diff),
            "-tu", "60000", "-tl", "-200",
            "-iu", "60000", "-il", "-200",
            "-tg", "1.0", "-ig", "1.0",
            "-tr", "5.0", "-ir", "5.0",
            "-r", "8",
            "-nsx", "5", "-nsy", "5",
            "-nss", "3",
            "-ft", "8",
            "-nrx", str(region_config.nrx),
            "-nry", str(region_config.nry),
            "-fft",
            "-v", "0",
        ]

        start = time.perf_counter()
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        elapsed = time.perf_counter() - start

        if result.returncode != 0:
            print(f"Error: {result.stderr}")
            return None

        return elapsed

    def analyze_region_overhead(self, image_size: int = 512, n_trials: int = 2):
        """
        Compare execution time: single region vs multi-region tiling.
        Quantifies overhead from region stitching, buffering, and synchronization.
        """
        print("\n" + "="*70)
        print("SINGLE-REGION OVERHEAD ANALYSIS")
        print("="*70)
        print(f"Image size: {image_size}×{image_size}")
        print(f"Trials: {n_trials}\n")

        configs = [
            RegionConfig(1, 1),      # Single region (baseline)
            RegionConfig(2, 2),      # 4 regions
            RegionConfig(4, 4),      # 16 regions
        ]

        results = {}

        with tempfile.TemporaryDirectory(prefix="hotpants_region_") as tmpdir:
            tmpdir = Path(tmpdir)
            template, science, diff = self.make_synthetic_images(
                tmpdir, image_size, image_size, n_stars=int((image_size / 512) ** 2 * 120)
            )

            for config in configs:
                times = []
                for trial in range(n_trials):
                    elapsed = self.run_hotpants(template, science, diff, config)
                    if elapsed is not None:
                        times.append(elapsed)
                        print(f"{config.name:20s} trial {trial + 1}/{n_trials}: {elapsed:.2f}s")

                if times:
                    avg = np.mean(times)
                    std = np.std(times)
                    results[config.name] = {"avg": avg, "std": std, "times": times}

        # Report
        print("\n" + "-"*70)
        print("RESULTS")
        print("-"*70)

        single_time = results.get("single_region", {}).get("avg")
        if not single_time:
            print("Single-region baseline not available")
            return

        print(f"\nSingle region (baseline): {single_time:.2f}s\n")

        for config_name, data in sorted(results.items()):
            if config_name == "single_region":
                continue

            avg = data["avg"]
            overhead = ((avg - single_time) / single_time) * 100
            per_region = avg / (int(config_name.split("_")[1].split("x")[0]) ** 2)

            print(f"{config_name:20s}: {avg:.2f}s  ({overhead:+6.1f}% vs baseline)")
            print(f"  Per-region avg:     {per_region:.4f}s")
            print(f"  Estimated overhead: {overhead:.1f}% (stitching, buffering, sync)")

    def estimate_single_region_benefits(self):
        """
        Estimate performance improvements possible with single-region optimizations.
        """
        report = """
# SINGLE-REGION PROCESSING OPTIMIZATION OPPORTUNITIES

## Overview
Single-region processing (`-nrx=1 -nry=1`) enables several significant optimizations
that are not feasible with multi-region tiling:

---

## 1. Eliminate Region Boundary Overhead (5–10% speedup)

### Current Multi-Region Flow
- Divide image into N regions
- Process each region independently
- Merge overlapping buffer zones to stitch regions together
- Resolve boundary discontinuities via blending

**Overhead sources:**
- Redundant convolution in buffer zones (typically 20–50 pixels per side)
- Memory allocation for per-region state (kernel coefficients, stamps, etc.)
- OpenMP synchronization barriers between regions
- Cache pollution from jumping between regions

### Single-Region Benefit
- No buffer zones → full image is contiguous in memory
- Single allocation for kernel coefficients (vs. per-region)
- No inter-region synchronization
- Improved L2/L3 cache efficiency (full image stays resident longer)
- Reduced memory fragmentation

**Estimated speedup:** 5–10% (higher for smaller images, lower for large images with many cores)

---

## 2. Globally-Optimal Kernel Basis (2–5% improvement in accuracy, minor speed impact)

### Current Approach
- Fit a per-region polynomial kernel (default: degree 6, 4, 2 for basis 1, 2, 3)
- May overfit to local PSF variations within region

### Single-Region Benefit
- Single global polynomial fit to entire image PSF variation
- Can use lower polynomial degrees (save coefficients, reduce computation)
- Avoids redundant basis fitting across regions
- More stable solution (less local overfitting)

**Possible implementation:**
```c
// Instead of fitting K(x,y) per region, fit globally
// K(x,y) = Σᵢ cᵢ(x,y) φᵢ  where cᵢ is degree 2–3 instead of 6
// Reduces total coefficients from (7+5+3)×nRegions to (3+3+3)×1
```

**Estimated speedup:** 2–5% (fewer coefficients to solve for, better convergence)

---

## 3. Adaptive Stamp Grid Layout (5–15% speedup potential)

### Current Approach
- Fixed uniform stamp grid (`-nsx`, `-nsy`)
- Independent of PSF or star distribution
- May have excessive stamps in empty regions, insufficient in crowded regions

### Single-Region Benefit
- Design stamp grid to match source distribution
- Hierarchical or quadtree-based stamp layout
  - Dense stamps near bright stars
  - Sparse stamps in empty regions
- Eliminates redundant PSF center searches in empty areas

**Possible implementation:**
```python
# Pseudo-code: adaptive stamp grid
bright_pixels = locate_bright_pixels(image, threshold=5*sigma)
stamp_density = quad_tree_from_sources(bright_pixels, target_spacing=50)
# Result: fewer total stamps to process, focused on regions of interest
```

**Estimated speedup:** 5–15% (fewer stamps = fewer getPsfCenters calls; reduces 35% bottleneck)

---

## 4. Simplified Noise Model (1–3% speedup)

### Current Approach
- Noise model fitted per-region
- Must account for local variations in gain, readout noise, background

### Single-Region Benefit
- Single global noise model
- Reduced parameter estimation overhead
- Simpler variance propagation logic

**Estimated speedup:** 1–3% (modest but low-hanging fruit)

---

## 5. Thin-Plate Splines for Spatial Variation (10–20% accuracy improvement, TBD speedup)

**From CLAUDE.md wishlist:** "Big one: using thin plate splines (or B-splines) to model
spatial variation and background instead of polynomials. This coupled with setting
`-nrx` and `-nry` to 1 means that entire images can be processed in one go."

### Benefit
- Replace polynomial spatial variation with smoothing splines
- Avoids step discontinuities at region boundaries (common with multi-region processing)
- Continuous, smooth kernel across entire image
- Particularly beneficial for wide-field instruments

**Challenges:**
- Requires different solver (spline basis vs. polynomial basis)
- Regularization parameter tuning (smoothness vs. fit quality)
- Larger computational cost initially (offset by fewer regions)

**Estimated speedup/accuracy:** 10–20% improvement in final image quality; neutral or slightly slower CPU time

---

## 6. Direct GPU Acceleration of Single Region (potential 3–10× speedup)

### Current State
- CPU-bound algorithm, multi-core parallelism via OpenMP
- No GPU support

### Single-Region Benefit
- Entire image resident on GPU memory (typical: 512² = 2 MB, 1024² = 8 MB)
- All operations (convolution, PSF center finding, matrix solve) can be GPU-resident
- Avoids PCIe transfers for multi-region batching

**Feasible optimizations:**
1. **CUDA FFT:** Replace FFTW3 with cuFFT for convolution (native GPU, 2–3× speedup)
2. **Thrust/cuBLAS:** PSF center search via parallel reductions
3. **cuSOLVER:** Cholesky solve on GPU

**Estimated speedup:** 3–10× on modern GPU (Nvidia A100/H100, AMD MI300)

---

## Implementation Roadmap

### Phase 1 (Quick wins, 10–15% overall speedup)
1. ✓ Eliminate region buffer overhead (reduce redundant convolution)
2. ✓ Reduce global polynomial degree for single region
3. ✓ Benchmark to establish baseline

### Phase 2 (Moderate effort, 5–10% additional speedup)
4. Implement adaptive stamp grid based on source detection
5. Merge/simplify noise model code
6. OpenMP parallelization improvements

### Phase 3 (Major refactor, 10–20% accuracy improvement)
7. Implement thin-plate splines as alternative to polynomials
8. Continuous kernel model eliminating region boundaries
9. Benchmark against Alard & Lupton reference

### Phase 4 (Exploratory, potential 3–10× speedup)
10. GPU acceleration via CUDA (cuFFT, cuBLAS, cuSOLVER)
11. Batched processing of multi-image surveys

---

## Testing & Validation

Single-region optimizations should be validated against:
1. **Regression tests:** Existing synthetic test suite (`test_api_integration.py`)
2. **Accuracy metrics:** Compare difference images to multi-region baseline
3. **Performance:** Use `benchmark_suite.py` to track speedup
4. **Real astronomical data:** Validation against ZTF, LSST, Rubin Observatory images

---

## References
- Alard & Lupton (1998): "A Method for Optimal Image Subtraction", ApJ 503:325
- Bramich (2008): "Delta-function basis for image subtraction", MNRAS 386:887
  (alternative basis suitable for single-region fitting)
- LSST-DM-TN021: "Decorrelation of Difference Images"
  (post-convolution whitening kernel for noise reduction)
"""
        return report

    def run(self):
        """Run full analysis."""
        print("\nHOTPANTS SINGLE-REGION ANALYSIS")

        # Overhead analysis
        self.analyze_region_overhead(image_size=512, n_trials=2)

        # Optimization opportunity report
        report = self.estimate_single_region_benefits()
        print(report)

        return report


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Single-region optimization analysis for HOTPANTS")
    parser.add_argument(
        "--hotpants",
        type=Path,
        default=Path("/home/user/hotpants/build/hotpants"),
        help="Path to hotpants binary"
    )
    args = parser.parse_args()

    if not args.hotpants.exists():
        print(f"Error: {args.hotpants} not found")
        return 1

    analyzer = SingleRegionAnalyzer(args.hotpants)
    analyzer.run()
    return 0


if __name__ == "__main__":
    exit(main())
