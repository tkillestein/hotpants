# Migration Guide: Global State → Context-Based API

As of May 2026, HOTPANTS introduces a **context-based API** that eliminates global state, improving thread safety and code clarity. This guide documents the transition.

## Background

Previously, HOTPANTS Python API used a **global_state() context manager** to temporarily set C library globals before calling functions. This approach worked but had limitations:

- **Thread safety**: Multiple threads accessing globals simultaneously could interfere
- **Transparency**: Hidden global modifications made debugging harder
- **Explicitness**: Resource lifecycle was implicit in the context manager

The new **context-based API** passes configuration and state explicitly through function calls, following a cleaner design pattern.

---

## Old API (Deprecated but Still Supported)

```python
from hotpants import fit_kernel, spatial_convolve, KernelConfig, RegionLayout, NoiseThresholds

# Load images
template = np.random.randn(256, 256).astype(np.float32)
science = np.random.randn(256, 256).astype(np.float32)

# Configure
config = KernelConfig(kernel_half_width=10, kernel_order=2)
layout = RegionLayout(n_regions_x=1, n_regions_y=1, stamps_per_region_x=10, stamps_per_region_y=10)
thresholds = NoiseThresholds(template_gain=1.0, science_gain=1.0)

# Use global_state() context manager (OLD WAY)
from hotpants._core import global_state

with global_state({'hwKernel': 10, 'kerOrder': 2, ...}):
    solution = fit_kernel(template, science, config=config, layout=layout, thresholds=thresholds)
    diff = spatial_convolve(science, solution, config)
```

**Note:** This still works but is deprecated. It may be removed in a future major version.

---

## New API (Recommended)

### Option 1: Use High-Level Context-Based Function

```python
from hotpants import fit_kernel_with_context, KernelConfig, RegionLayout, NoiseThresholds
import numpy as np

template = np.random.randn(256, 256).astype(np.float32)
science = np.random.randn(256, 256).astype(np.float32)

config = KernelConfig(kernel_half_width=10, kernel_order=2)
layout = RegionLayout(n_regions_x=1, n_regions_y=1)
thresholds = NoiseThresholds(template_gain=1.0, science_gain=1.0)

# Call context-based function (NEW WAY - recommended)
solution = fit_kernel_with_context(
    template, science,
    config=config,
    layout=layout,
    thresholds=thresholds,
    verbose=1,
    n_thread=4
)

print(f"Kernel fit chi2: {solution.chi2:.3f}")
print(f"Kernel norm: {solution.kernel_norm:.3f}")
```

**Benefits:**
- No global state manipulation
- Thread-safe: each call has its own context
- Clear resource ownership and cleanup
- Explicit configuration passing

### Option 2: Manual Context Management (Advanced)

For fine-grained control, you can manage contexts explicitly:

```python
from hotpants._core import (
    create_config_struct,
    create_kernel_context,
    create_region_context,
    cleanup_kernel_context,
    cleanup_region_context,
    fitkernel_v2,
    spatial_convolve_v2,
    allocate_stamps,
    free_stamps,
)

template = np.random.randn(256, 256).astype(np.float32)
science = np.random.randn(256, 256).astype(np.float32)

# Create configuration struct
config = create_config_struct()
config.kernel_half_width = 10
config.kernel_order = 2
config.n_regions_x = 1
config.n_regions_y = 1
config.stamps_per_region_x = 10
config.stamps_per_region_y = 10

# Create kernel context (allocates basis vectors once)
ny, nx = template.shape
kctx = create_kernel_context(
    config, nx, ny,
    config.n_regions_x, config.n_regions_y,
    config.stamps_per_region_x, config.stamps_per_region_y
)

# Create region context (allocates per-region buffers)
n_stamps = config.stamps_per_region_x * config.stamps_per_region_y
rctx = create_region_context(kctx, config, n_stamps)

try:
    # Allocate stamps
    stamps = allocate_stamps(n_stamps)

    try:
        # Call algorithm with explicit contexts
        kernel_coeffs, meansig, scatter, n_skipped = fitkernel_v2(
            stamps, n_stamps,
            template, science, None,  # no noise image
            config, kctx, rctx
        )

        # Convolve with fitted kernel
        diff_image, conv_mask = spatial_convolve_v2(
            science, kernel_coeffs,
            config, kctx, rctx
        )

    finally:
        free_stamps(stamps, n_stamps)

finally:
    cleanup_region_context(rctx)
    cleanup_kernel_context(kctx)
```

**This approach is useful for:**
- Reusing contexts across multiple operations (efficiency)
- Custom control flow and error handling
- Integration with your own multi-threading code

---

## Migrating Existing Code

### Step 1: Replace global_state() with fit_kernel_with_context()

**Before:**
```python
with global_state({'hwKernel': 10, 'kerOrder': 2}):
    solution = fit_kernel(template, science)
```

**After:**
```python
solution = fit_kernel_with_context(template, science, config=config)
```

### Step 2: Update Configuration

The new API uses **Pydantic dataclasses** (KernelConfig, RegionLayout, NoiseThresholds) instead of passing globals.

**Before:**
```python
with global_state({
    'hwKernel': 10,
    'kerOrder': 2,
    'bgOrder': 1,
    'tGain': 1.0,
    'iGain': 1.0,
    ...
}):
    solution = fit_kernel(template, science)
```

**After:**
```python
config = KernelConfig(kernel_half_width=10, kernel_order=2, bg_order=1)
thresholds = NoiseThresholds(template_gain=1.0, science_gain=1.0)
solution = fit_kernel_with_context(
    template, science,
    config=config,
    thresholds=thresholds
)
```

### Step 3: Thread-Safe Multi-Image Processing

**Old approach (NOT thread-safe):**
```python
# DON'T DO THIS IN PARALLEL — global state will interfere
for image_pair in image_pairs:
    with global_state({'hwKernel': 10}):
        solution = fit_kernel(...)
```

**New approach (thread-safe):**
```python
from concurrent.futures import ThreadPoolExecutor

def process_image_pair(template, science):
    return fit_kernel_with_context(template, science)

with ThreadPoolExecutor(max_workers=4) as executor:
    solutions = executor.map(process_image_pair, templates, sciences)
```

Each thread now has its own context — no interference.

---

## Compatibility

| Feature | Old API (global_state) | New API (context) |
|---------|------------------------|-------------------|
| **Thread safety** | ❌ No | ✅ Yes |
| **Global modifications** | ✅ Yes | ❌ No |
| **Explicit resource lifecycle** | ❌ Hidden | ✅ Clear |
| **Performance** | ✅ Same | ✅ Same |
| **Backwards compatible** | — | ✅ Old API still works |

---

## FAQ

### Q: Will my old code break?

**A:** No. The old `fit_kernel()` and `global_state()` are still available. Deprecation warnings will appear in future versions, but code will continue to work.

### Q: Is there a performance difference?

**A:** No. The new context API has the same performance. Context initialization is negligible compared to kernel fitting time.

### Q: Can I mix old and new APIs?

**A:** Not recommended. Mixing globals (old API) with explicit contexts (new API) may cause unexpected behavior. Choose one approach per operation.

### Q: How do I use this in a package that calls HOTPANTS?

**A:** Wrap the call in `fit_kernel_with_context()`:

```python
def my_image_differencing_function(template, science, config=None):
    if config is None:
        config = KernelConfig()
    return fit_kernel_with_context(template, science, config=config)
```

### Q: What if I need more control?

**A:** Use manual context management (Option 2 above) for fine-grained control over resource lifecycle and multi-threaded scenarios.

---

## References

- **Context Structures:** `src/hotpants_context.h`
- **Python Wrappers:** `src/hotpants/_core.py`
- **C Implementation:** `src/hotpants_wrapper.c`, `src/hotpants_context.c`
- **Architecture:** See CLAUDE.md section "Global State Elimination" for design details

---

**Questions?** See CLAUDE.md or file an issue at https://github.com/tkillestein/hotpants/issues.
