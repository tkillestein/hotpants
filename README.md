# HOTPANTS

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Build Status](https://img.shields.io/github/actions/workflow/status/tkillestein/hotpants/ci.yml?branch=master)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)
![C17](https://img.shields.io/badge/C-17-blue.svg)

**High Order Transform of PSF ANd Template Subtraction**

A fast, production-grade tool for image differencing. HOTPANTS fits a spatially-varying
convolution kernel to match two astronomical images' point-spread functions (PSFs),
producing optimal difference images for transient detection.

> Based on the algorithm
> from [Alard & Lupton (1998)](https://iopscience.iop.org/article/10.1086/305984).
> HOTPANTS was written by Andy Becker; maintained and modernized here.

---

## Documentation

📖 **[Full Documentation](https://tkillestein.github.io/hotpants/)** on GitHub Pages

- API reference (C)
- Installation and quick-start guides
- PSF matching strategy and tuning
- Algorithm overview and references

---

## Features

- **Fast convolution** — FFT-accelerated via FFTW3 (3–8× speedup on typical images)
- **Parallelized** — multi-threaded with OpenMP for modern CPUs (4–8× total speedup)
- **Decorrelation** — LSST DMTN-021 afterburner noise decorrelation with spatially-varying whitening kernel
- **Modern build system** — CMake with automatic dependency detection
- **Comprehensive documentation** — Doxygen-generated API docs and user guides
- **Python API** — full numpy integration via ctypes for the scientific Python stack

---

## Installation

### From Source (Recommended for Development)

1. **Install system dependencies:**

   **macOS:**
   ```bash
   brew install cfitsio openblas fftw lapack
   ```

   **Ubuntu/Debian:**
   ```bash
   sudo apt-get install libcfitsio-dev libopenblas-dev libfftw3-dev liblapacke-dev
   ```

2. **Clone and build C library:**
   ```bash
   git clone https://github.com/tkillestein/hotpants.git
   cd hotpants
   cmake -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j$(nproc)
   ```

3. **Install Python package:**
   ```bash
   pip install -e .
   # Or with development dependencies:
   pip install -e ".[dev]"
   ```

4. **Set library path (if not in system path):**
   ```bash
   export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH
   ```

See [CONTRIBUTING.md](CONTRIBUTING.md) for full development setup, testing, and profiling instructions.

### Build Configuration

The CMake build auto-detects CFITSIO, OpenBLAS, FFTW3, and OpenMP:
- FFTW3 is **required** (no fallback to direct convolution)
- OpenMP is optional; builds without it if unavailable

Optional CMake flags:
- `-DUSE_OPENMP=ON/OFF` — enable/disable multi-threading (default: auto-detect)
- `-DCMAKE_C_FLAGS="-O3 -march=native"` — optimize for local machine

### Troubleshooting

- **"Could not find libhotpants":** Set `LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH`
- **CMake fails:** Ensure all dependencies are installed (CFITSIO, OpenBLAS, FFTW3, LAPACK)
- **Python import errors:** C library must be built first (see Installation)

---

## Quick Start

```bash
# Basic usage: match template to science image
hotpants -inim science.fits -tmplim template.fits -outim difference.fits

# Use 4 threads
OMP_NUM_THREADS=4 hotpants -inim science.fits -tmplim template.fits -outim difference.fits

# Parallel processing over image regions
hotpants -inim science.fits -tmplim template.fits -outim difference.fits \
  -nrx 2 -nry 2

# Save kernel info, output noise map
hotpants -inim science.fits -tmplim template.fits -outim difference.fits \
  -savexy kernel_positions.txt -oni noise.fits
```

### PSF Matching Strategy

The choice of convolution direction depends on your PSF sizes:

- **Science PSF sharper than template**: Convolution can introduce artifacts. Options:
    - Pre-convolve science image with its PSF before matching (recommended)
    - Force convolution on template: `-c t`

- **Science PSF broader than template** (common): Match succeeds naturally. Set basis
  Gaussians centered on `σ = √(σ²_science − σ²_template)`:
  ```bash
  hotpants ... -ng 3 6 0.5*sigma 4 sigma 2 2.0*sigma
  ```

See [CLAUDE.md](CLAUDE.md) for algorithm details and performance tuning.

### Python API

```python
import numpy as np
from hotpants import fit_kernel, spatial_convolve, KernelConfig

# Load images (float32)
template = np.load('template.npy').astype(np.float32)
science = np.load('science.npy').astype(np.float32)

# Configure and fit kernel
config = KernelConfig(kernel_half_width=15, kernel_order=2)
kernel_solution = fit_kernel(template, science, config=config)

# Create difference image
difference = spatial_convolve(science, kernel_solution)
np.save('difference.npy', difference)
```

**Key classes:**
- `KernelConfig` — kernel fitting parameters
- `RegionLayout` — image tiling configuration
- `NoiseThresholds` — data quality thresholds
- `KernelSolution` — fitted kernel coefficients and metrics

**Requirements:** C library must be built first; numpy arrays must be float32 and same shape.

---

## Output

HOTPANTS produces a multi-layer FITS file:

| Layer | Content                                                                   |
|-------|---------------------------------------------------------------------------|
| 1     | Difference image: `D = I − T⊗K`                                           |
| 2     | Decorrelated difference image: `D_decorr = φ ⊗ D` (if `-useDecorrelation 1`) |
| 3+    | Optional: convolved image, noise map, noise-scaled difference, pixel mask |

Use `-allm` to output all available layers. Use `-nim`, `-oni` for individual layers.

### Decorrelation Output

When decorrelation is enabled (`-useDecorrelation 1`), the output FITS file includes:

- **Layer 2 (DECORR):** Decorrelated difference image with correlated noise removed
- **Header keywords:**
  - `HIERARCH DECORR_METHOD`: "AFTERBURNER" (LSST DMTN-021)
  - `HIERARCH DECORR_SCI_VAR`: Science image variance used
  - `HIERARCH DECORR_TPL_VAR`: Template image variance used

The decorrelated image is obtained by convolving D with a spatially-varying whitening kernel
φ(x,y) that removes the correlation structure introduced by the PSF-matching convolution.

---

## Configuration

Full option reference: run `hotpants` with no arguments.

Key tuning parameters:

| Option         | Purpose                                     | Default |
|----------------|---------------------------------------------|---------|
| `-r`           | Convolution kernel half-width (pixels)      | 10      |
| `-nrx`, `-nry` | Region grid (for spatial kernel variation)  | 1×1     |
| `-nsx`, `-nsy` | Stamp grid within each region               | 10×10   |
| `-ko`          | Spatial polynomial order (kernel variation) | 2       |
| `-bgo`         | Spatial polynomial order (background)       | 1       |
| `-ft`          | Centroid fitting threshold (RMS pixels)     | 20.0    |
| `-ks`          | Bad-stamp rejection threshold (sigma)       | 2.0     |

### Decorrelation Parameters

| Option                   | Purpose                                | Default |
|--------------------------|----------------------------------------|---------|
| `-useDecorrelation`      | Enable DMTN-021 decorrelation          | 0       |
| `-decorrScienceVar`      | Science image variance (auto-estimated if 0) | 0 |
| `-decorrTemplateVar`     | Template image variance (auto-estimated if 0) | 0 |
| `-decorrUseTPS`          | Use TPS interpolation for φ field      | 0       |

**When to use decorrelation:**
- Image surveys with time-domain analysis (transient detection, variability)
- When correlated noise in difference images affects photometry
- Production pipelines requiring optimal signal-to-noise on faint transients
- Single-region mode (`-nrx 1 -nry 1`) for best results

---

## Performance

On a modern 4-core CPU with FFTW acceleration:

| Image size | Time  | Speedup vs. direct |
|------------|-------|--------------------|
| 1024×1024  | ~0.1s | 3–4×               |
| 4096×4096  | ~2s   | 5–8×               |

Use `OMP_NUM_THREADS=N` to control threading. Typical scaling: 2–4× for 4 threads.

---

## Building from Source

### Requirements

- **CFITSIO** ≥ 3.45 (FITS I/O)
- **BLAS/LAPACK** (linear algebra)
- **FFTW3** ≥ 3.3 (mandatory for FFT convolution)
- OpenMP (optional, for multi-threading)
- CMake ≥ 3.10

### Compiler Optimization

For maximum performance:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-O3 -march=native -funroll-loops"
```

Use `-march=native` for local builds; omit for portable binaries.

---

## Development

See [CLAUDE.md](CLAUDE.md) for:

- Algorithm details and code structure
- Completed optimizations (FFT, SIMD vectorization, cache-aware tiling, parallelism)
- Performance profiling workflow
- Technical debt and refactoring opportunities

See [CONTRIBUTING.md](CONTRIBUTING.md) for build, testing, and profiling guidelines.

---

## License

MIT License, 2013–2026. See [LICENSE](LICENSE).

## Citation

If you use this fork of HOTPANTS in research, please cite the original algorithm and the
original
HOTPANTS implementation:

```bibtex
@article{Alard1998,
  author = {Alard, C. and Lupton, R. H.},
  title = {A Method for Optimal Image Subtraction},
  journal = {ApJ},
  volume = {503},
  pages = {325},
  year = {1998},
  doi = {10.1086/305984}
}

@misc{Becker2015,
  author = {Becker, A. C.},
  title = {HOTPANTS: High Order Transform of PSF ANd Template Subtraction},
  year = {2015},
  url = {https://ascl.net/1504.004},
  note = {Astrophysics Source Code Library, record ascl:1504.004}
}
```

---

**Questions?** Open an issue on [GitHub](https://github.com/tkillestein/hotpants).
