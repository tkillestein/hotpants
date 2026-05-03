"""
HOTPANTS — High Order Transform of PSF ANd Template Subtraction

A Python API for fitting and applying spatially-varying convolution kernels
for astronomical image differencing.

Core Reference: Alard & Lupton (1998), ApJ 503:325
https://iopscience.iop.org/article/10.1086/305984

Example usage:

    import numpy as np
    from hotpants import fit_kernel, spatial_convolve, KernelConfig, RegionLayout, NoiseThresholds

    # Load images (assumed numpy arrays, float32)
    template = np.load('template.npy')
    science = np.load('science.npy')

    # Configure kernel fitting
    config = KernelConfig()
    layout = RegionLayout()
    thresholds = NoiseThresholds()

    # Fit spatially-varying kernel
    kernel_solution = fit_kernel(template, science, None, config, layout, thresholds)

    # Apply kernel to science image to create difference image
    difference = spatial_convolve(science, kernel_solution, config)
"""

__version__ = "0.0.1"
__author__ = "Andy Becker (C core), Python bindings"
__license__ = "MIT"

# Deferred import of API functions to avoid cffi build errors during package discovery
def __getattr__(name):
    """Lazy load public API on first access."""
    if name in ("fit_kernel", "spatial_convolve", "KernelConfig",
                "RegionLayout", "NoiseThresholds", "KernelSolution"):
        from . import api
        return getattr(api, name)
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

__all__ = [
    "fit_kernel",
    "spatial_convolve",
    "KernelConfig",
    "RegionLayout",
    "NoiseThresholds",
    "KernelSolution",
]
