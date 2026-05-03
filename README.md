# HOTPANTS

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Build Status](https://img.shields.io/github/actions/workflow/status/tkillestein/hotpants/ci.yml?branch=main)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)
![C99](https://img.shields.io/badge/C-99-blue.svg)

**High Order Transform of PSF ANd Template Subtraction**

A fast, production-grade tool for image differencing. HOTPANTS fits a spatially-varying convolution kernel to match two astronomical images' point-spread functions (PSFs), producing optimal difference images for transient detection.

> Based on the algorithm from [Alard & Lupton (1998)](https://iopscience.iop.org/article/10.1086/305984). Originally written by Andy Becker; maintained and modernized here.

---

## Documentation

📖 **[Full Documentation](https://tkillestein.github.io/hotpants/)** on GitHub Pages
- API reference (C)
- Installation and quick-start guides
- PSF matching strategy and tuning
- Algorithm overview and references

---

## Features

- **Fast convolution** — FFT-accelerated via FFTW3 for large images
- **Parallelized** — multi-threaded with OpenMP for modern CPUs
- **Flexible kernel** — spatially-varying polynomial basis, customizable basis functions
- **Robust masking** — propagates bad pixel masks through the pipeline
- **Production-ready** — ~60% reduction in difference image noise vs. simple differencing

---

## Installation

### CMake (recommended)

```bash
# Install dependencies (macOS)
brew install cfitsio openblas fftw

# Install dependencies (Ubuntu/Debian)
sudo apt-get install libcfitsio-dev libopenblas-dev libfftw3-dev

# Build
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
sudo cmake --install build  # optional
```

The CMake build auto-detects CFITSIO, OpenBLAS, FFTW3, and OpenMP. Use `-DUSE_FFTW=OFF` to disable FFT acceleration if needed.

---

## Quick Start

```bash
# Basic usage: match template to science image
hotpants -inim science.fits -tmplim template.fits -outim difference.fits

# Enable FFT-accelerated convolution (much faster for large images)
hotpants -inim science.fits -tmplim template.fits -outim difference.fits -fft

# Use 4 threads
OMP_NUM_THREADS=4 hotpants -inim science.fits -tmplim template.fits -outim difference.fits

# Parallel processing over image regions
hotpants -inim science.fits -tmplim template.fits -outim difference.fits \
  -nrx 2 -nry 2 -fft

# Save kernel info, output noise map
hotpants -inim science.fits -tmplim template.fits -outim difference.fits \
  -savexy kernel_positions.txt -oni noise.fits
```

### PSF Matching Strategy

The choice of convolution direction depends on your PSF sizes:

- **Science PSF sharper than template**: Convolution can introduce artifacts. Options:
  - Pre-convolve science image with its PSF before matching (recommended)
  - Force convolution on template: `-c t`
  
- **Science PSF broader than template** (common): Match succeeds naturally. Set basis Gaussians centered on `σ = √(σ²_science − σ²_template)`:
  ```bash
  hotpants ... -ng 3 6 0.5*sigma 4 sigma 2 2.0*sigma
  ```

See [CLAUDE.md](CLAUDE.md) for algorithm details and performance tuning.

---

## Output

HOTPANTS produces a multi-layer FITS file:

| Layer | Content |
|-------|---------|
| 1 | Difference image: `D = I − T⊗K` |
| 2+ | Optional: convolved image, noise map, noise-scaled difference, pixel mask |

Use `-allm` to output all available layers. Use `-nim`, `-cim`, `-oni` for individual layers.

---

## Configuration

Full option reference: run `hotpants` with no arguments.

Key tuning parameters:

| Option | Purpose | Default |
|--------|---------|---------|
| `-r` | Convolution kernel half-width (pixels) | 10 |
| `-nrx`, `-nry` | Region grid (for spatial kernel variation) | 1×1 |
| `-nsx`, `-nsy` | Stamp grid within each region | 10×10 |
| `-ko` | Spatial polynomial order (kernel variation) | 2 |
| `-bgo` | Spatial polynomial order (background) | 1 |
| `-ft` | Centroid fitting threshold (RMS pixels) | 20.0 |
| `-ks` | Bad-stamp rejection threshold (sigma) | 2.0 |
| `-fft` | Enable FFT acceleration | off (use direct convolution) |

---

## Performance

On a modern 4-core CPU with FFTW acceleration:

| Image size | Time | Speedup vs. direct |
|------------|------|-------------------|
| 1024×1024 | ~0.1s | 3–4× |
| 4096×4096 | ~2s | 5–8× |

Use `OMP_NUM_THREADS=N` to control threading. Typical scaling: 2–4× for 4 threads.

---

## Building from Source

### Requirements

- CFITSIO ≥ 3.45 (FITS I/O)
- OpenBLAS or LAPACK (linear algebra)
- FFTW3 ≥ 3.3 (optional, for FFT acceleration)
- OpenMP (optional, for multi-threading)
- CMake ≥ 3.10 or GNU Make

### Compiler Optimization

For maximum performance on your machine:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_FLAGS="-O3 -march=native -funroll-loops"
```

Use `-march=native` for local builds; omit for portable binaries.

---

## Development

See [CLAUDE.md](CLAUDE.md) for:
- Algorithm overview and modernization roadmap
- Build system details (CMake)
- Performance profiling
- Planned optimizations (SIMD, adaptive basis, GPU offload)

Contributions welcome! See future `CONTRIBUTING.md` for guidelines.

---

## License

MIT License, 2013–2026. See [LICENSE](LICENSE).

## Citation

If you use HOTPANTS in research, please cite the original algorithm and this implementation:

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
