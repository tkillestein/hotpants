# Contributing to HOTPANTS

Thank you for interest in contributing to HOTPANTS! This guide covers development setup, building, testing, and submitting changes.

## Development Setup

### Requirements

- **Python:** ≥ 3.11
- **CMake:** ≥ 3.20
- **C Compiler:** gcc/clang with C17 support
- **System Libraries:**
  - CFITSIO (`libcfitsio-dev` on Ubuntu/Debian, `cfitsio` on macOS)
  - OpenBLAS (`libopenblas-dev` on Ubuntu/Debian, `openblas` on macOS)
  - FFTW3 (`libfftw3-dev` on Ubuntu/Debian, `fftw` on macOS)
  - LAPACK with LAPACKE (`liblapacke-dev` on Ubuntu/Debian, included in Accelerate on macOS)

### Installation

**macOS:**
```bash
brew install cfitsio openblas fftw lapack
```

**Ubuntu/Debian:**
```bash
sudo apt-get install libcfitsio-dev libopenblas-dev libfftw3-dev liblapacke-dev
```

**Then:**
```bash
git clone https://github.com/tkillestein/hotpants.git
cd hotpants
uv sync --dev  # Install Python dependencies
```

## Building the C Library

The C library must be built before running tests or using the Python API:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

This creates `build/libhotpants.so` (Linux) or `build/libhotpants.dylib` (macOS).

**Optional CMake flags:**
- `-DUSE_FFTW=ON/OFF` — Enable/disable FFT convolution (default: auto-detect)
- `-DUSE_OPENMP=ON/OFF` — Enable/disable OpenMP parallelism (default: auto-detect)
- `-DCMAKE_BUILD_TYPE=Debug/Release` — Build type (default: Release)

## Running Tests

**All tests:**
```bash
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH
uv run pytest tests/ -v
```

**Specific test file:**
```bash
uv run pytest tests/test_api.py -v
```

**Single test:**
```bash
uv run pytest tests/test_api.py::TestKernelConfig::test_default_config -v
```

**With coverage:**
```bash
uv run pytest tests/ --cov=src/hotpants --cov-report=html
```

## Code Style & Conventions

### Python Code

- **Style:** PEP 8, enforced by Ruff
- **Type hints:** All public functions must have type annotations
- **Docstrings:** Google-style docstrings for public APIs
- **Line length:** 100 characters (configured in pyproject.toml)

**Check style:**
```bash
uv run ruff check src/hotpants tests
```

**Auto-fix style issues:**
```bash
uv run ruff check --fix src/hotpants tests
```

### C Code

- **Style:** Modern C17 (see `CMakeLists.txt`)
- **Naming:** snake_case for functions, UPPERCASE for constants
- **Documentation:** Doxygen comments for public API (`@brief`, `@param`, `@return`)
- **Global variables:** Follow prefix convention in globals.h (t*, i*, m*, o*)
- **Magic numbers:** Use #define constants in defaults.h

See [CLAUDE.md](CLAUDE.md) for algorithm references and implementation details.

## Making Changes

### 1. Create a feature branch

```bash
git checkout -b feature/your-feature-name
```

### 2. Make changes and commit

```bash
git add src/
git commit -m "Brief description of changes

More detailed explanation if needed. Reference related issues or papers.

https://claude.ai/code/session_01CZEfmSpumdfwufTRmH6C1y"
```

**Commit message guidelines:**
- First line: ≤50 characters, imperative mood ("Add feature", not "Added feature")
- Blank line, then detailed explanation if needed
- Reference issues/PRs: "Fixes #123" or "Related to #456"
- Include session link for AI-assisted work

### 3. Run tests locally

```bash
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH
uv run pytest tests/ -v
```

### 4. Push and create a PR

```bash
git push -u origin feature/your-feature-name
```

Then open a pull request on GitHub. Include:
- Description of changes
- Test results
- Any performance impact
- Relevant issue references

## Performance Profiling

### Baseline (single-threaded)
```bash
OMP_NUM_THREADS=1 ./build/hotpants -inim science.fits -tmplim template.fits -outim diff.fits
```

### Profile with gprof
```bash
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build

OMP_NUM_THREADS=1 ./build/hotpants [options]
gprof ./build/hotpants gmon.out | head -50
```

### Profile with perf
```bash
perf record -g ./build/hotpants [options]
perf report
```

### Flame graphs
```bash
perf record -F 99 -g ./build/hotpants [options]
perf script > out.perf
# Install flamegraph: git clone https://github.com/brendangregg/FlameGraph
flamegraph.pl out.perf > graph.svg
```

## Documentation

### Building Sphinx docs
```bash
cd docs
make html
open _build/html/index.html
```

### Adding to documentation

- **API docs:** Edit C code Doxygen comments in `src/`
- **User guides:** Add `.rst` files to `docs/guides/`
- **Examples:** Add `.py` files to `docs/examples/`

After editing, rebuild and commit the changes.

## Troubleshooting

### "Could not find libhotpants"

The C library wasn't found. Try:
```bash
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH
```

Or build the library:
```bash
cmake -B build && cmake --build build
```

### CMake configuration fails

Missing system dependencies. Run the appropriate install command:
- **macOS:** `brew install cfitsio openblas fftw`
- **Ubuntu:** `sudo apt-get install libcfitsio-dev libopenblas-dev libfftw3-dev liblapacke-dev`

### Tests timeout or hang

Check if OpenMP is causing issues:
```bash
OMP_NUM_THREADS=1 uv run pytest tests/ -v
```

## Priorities for Contributors

From [CLAUDE.md](CLAUDE.md):

1. **High priority:** Python integration tests (Phase 3)
2. **Medium priority:** Phase 2 readability improvements (loop variable refactoring, function decomposition)
3. **Future work:** Bramich kernel basis, whitening kernels, thin-plate spline spatial modeling

See [CLAUDE.md § Improvements Wishlist](CLAUDE.md#improvements-wishlist) for more ideas.

## Questions?

- Check [CLAUDE.md](CLAUDE.md) for algorithm and architecture details
- Review existing PRs for examples of accepted contributions
- Open an issue to discuss larger changes before implementing

Thank you for contributing!
