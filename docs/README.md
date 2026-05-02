# Building Documentation

HOTPANTS documentation is built with **Sphinx** and **Doxygen**.

## Prerequisites

- `sphinx >= 7.0`
- `breathe >= 4.35`
- `sphinx-rtd-theme >= 2.0`
- `doxygen >= 1.9` (for C API documentation)

## Installation

```bash
# Install dependencies via uv
uv pip install sphinx breathe sphinx-rtd-theme

# Install Doxygen (system package)
# macOS
brew install doxygen

# Linux (Ubuntu/Debian)
sudo apt-get install doxygen
```

## Building Locally

```bash
# From the docs directory
cd docs

# Build HTML documentation (includes Doxygen + Sphinx)
make html

# View in browser
open _build/html/index.html  # macOS
xdg-open _build/html/index.html  # Linux

# Clean build artifacts
make clean
```

## Documentation Structure

```
docs/
├── conf.py                # Sphinx configuration
├── index.rst              # Main documentation index
├── Doxyfile               # Doxygen configuration (C API)
├── guides/
│   ├── installation.rst
│   ├── quickstart.rst
│   ├── psf_matching.rst
│   ├── tuning.rst
│   └── algorithm.rst
├── api/
│   ├── c_api.rst          # C API (auto-generated from Doxygen)
│   └── python_api.rst     # Python API (placeholder for future bindings)
└── _build/                # Build artifacts (git-ignored)
```

## Continuous Integration

Documentation is automatically built and deployed to GitHub Pages on every push to `main`. See `.github/workflows/docs.yml`.

## Workflow

1. **Write/edit documentation** in `.rst` (reStructuredText) or Markdown
2. **Add Doxygen comments** to C functions in `alard.c` and `functions.c`
3. **Build locally** to verify (requires Doxygen + Sphinx)
4. **Push to GitHub**; CI/CD automatically rebuilds and deploys

## Doxygen Comments

Doxygen comments use the JavaDoc-style syntax:

```c
/**
 * @brief One-line summary.
 *
 * @details Longer description, can span multiple paragraphs and reference
 * equations or papers (e.g., Alard & Lupton 1998, Eq. 2).
 *
 * @param[in] param_in Input parameter description.
 * @param[out] param_out Output parameter description.
 * @param[in,out] param_both Input/output parameter.
 * @return Description of return value.
 *
 * @note Any special notes or invariants.
 * @see Related function or section.
 */
```

See existing comments in `alard.c` and `functions.c` for examples.

## Troubleshooting

**"WARNING: Could not find reference target"**

Usually means a function isn't documented by Doxygen (not extracted). Check:
- Function is not static (unless EXTRACT_STATIC is set in Doxyfile)
- Doxygen comment is properly formatted
- Function signature is in `functions.h` or similar

**"ERROR: Breaker format error"**

Usually means malformed reStructuredText in `.rst` files. Check:
- Underline characters match (= or -) header length
- Code blocks properly indented
- No stray special characters

**Sphinx build fails**

1. Verify Sphinx is installed: `sphinx-build --version`
2. Check `conf.py` syntax: `python -m py_compile docs/conf.py`
3. Rebuild from scratch: `make clean && make html`

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md) (coming soon) for documentation guidelines.

## Resources

- [Sphinx Documentation](https://www.sphinx-doc.org/)
- [Breathe Documentation](https://breathe.readthedocs.io/)
- [Doxygen Manual](https://doxygen.nl/)
- [reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
