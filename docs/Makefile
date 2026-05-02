# Sphinx documentation build Makefile

.PHONY: help clean html doxygen

SPHINXBUILD   = sphinx-build
SPHINXOPTS    = -W --keep-going
SOURCEDIR     = .
BUILDDIR      = _build

help:
	@echo "Targets:"
	@echo "  doxygen  - Generate Doxygen XML documentation"
	@echo "  html     - Build HTML documentation with Sphinx"
	@echo "  clean    - Remove build artifacts"

doxygen:
	@echo "Building Doxygen XML..."
	doxygen ../Doxyfile

html: doxygen
	@echo "Building Sphinx HTML..."
	$(SPHINXBUILD) -b html $(SPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/html
	@echo "HTML documentation built in $(BUILDDIR)/html/"

clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILDDIR) _doxygen
