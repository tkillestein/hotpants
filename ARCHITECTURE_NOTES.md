# Architecture Notes: Kernel Basis Abstraction

## Current State

HOTPANTS implements a pluggable kernel basis system (Gaussian vs Delta) via `basis_t` interface in `basis.h`, but the implementation has a **partial abstraction gap** in the stamp convolution loop.

## The Issue

### Interface Definition (basis.h)
The `kernel_basis_t` struct defines `convolve_stamp()` callback:
```c
int (*convolve_stamp)(double* stamp, int mStampX, int mStampY,
                      int basisIdx, double* pixFitted);
```

### Actual Implementation (fillStamp loop)

**Gaussian basis** — alard.c:1005-1020
```c
for (gaussianCompIdx = 0; gaussianCompIdx < ngauss; gaussianCompIdx++) {
    for (idegx = ...; idegx <= ...; idegx++) {
        for (idegy = ...; idegy <= ...; idegy++) {
            xy_conv_stamp(stamp, imConv, vectorComponentIdx, ...);  // ← DIRECT CALL
            ++vectorComponentIdx;
        }
    }
}
```

**Delta basis** — alard.c:1021-1033
```c
for (int basisIdx = 0; basisIdx < nCompKer; basisIdx++) {
    if (xy_conv_stamp_delta(stamp, imConv, basisIdx) < 0)  // ← DIRECT CALL
        return 1;
    ++vectorComponentIdx;
}
```

## Problem

1. **`convolve_stamp()` callback is NOT used** in either implementation
   - Gaussian: basis_gaussian.c marks it as "no-op, vectors pre-filled"
   - Delta: basis_delta.c implements the callback but fillStamp() doesn't call it

2. **Functions are hardcoded in fillStamp()** rather than dispatched through `active_basis`
   - Breaks encapsulation: fillStamp() knows about implementation details
   - Makes it hard to add new basis types (must modify fillStamp loop)

3. **Inconsistent patterns**
   - Gaussian and Delta use completely different approaches
   - Gaussian: pre-computed basis vectors, applied in convolution
   - Delta: basis functions are defined implicitly by kernel pixel indices

## Proper Generic Design Would Be

```c
void fillStamp(stamp_struct* stamp, float* imConv, float* imRef) {
    // ...
    
    /* Generic dispatch through basis interface */
    for (int basisIdx = 0; basisIdx < nCompKer; basisIdx++) {
        if (active_basis->convolve_stamp(
            stamp_pixel_data,           // input stamp
            fwKSStamp, fwKSStamp,       // dimensions
            basisIdx,                   // which basis function
            stamp->vectors[basisIdx]    // output vector
        ) < 0) {
            return 1;  // error
        }
    }
    
    /* ... rest of fillStamp (background polynomial) ... */
}
```

Then both `basis_gaussian_convolve_stamp()` and `basis_delta_convolve_stamp()` would be called through the uniform interface.

## Why It Wasn't Done This Way

1. **Gaussian basis predates delta basis abstraction**
   - Original code: xy_conv_stamp() was called directly
   - Delta basis added later with its own xy_conv_stamp_delta()
   - Abstraction layer (basis.h) created but not fully applied

2. **Gaussian basis has inherently different semantics**
   - Gaussian: Fixed set of basis functions (3 Gaussians + polynomial orders)
   - Application: Separable 2D convolution of precomputed filters
   - The filters are precomputed outside fillStamp() and reused
   
3. **Delta basis is simpler but different**
   - Delta: nCompKer = fwKernel² basis functions (one per pixel)
   - Application: Direct pixel extraction with shifted indices
   - No precomputation needed

## Recommendations

### Short Term (No Change Needed)
- Current implementation works correctly
- Basis dispatch logic (set_active_basis) is clean
- Only fillStamp() stamp convolution loop deviates from pattern
- Code is well-commented

### Long Term (Refactoring Opportunity)
If adding more basis types (PCA, Fourier, etc.):

1. **Refactor fillStamp() to use basis interface:**
   ```c
   /* After Gaussian precomputation, before stamp loop */
   for (int basisIdx = 0; basisIdx < nCompKer; basisIdx++) {
       active_basis->convolve_stamp(..., basisIdx, ...);
   }
   ```

2. **Gaussian basis would gain a no-op convolve_stamp:**
   - Mark with comment explaining it's pre-filled by xy_conv_stamp()
   - Or: Move Gaussian filter application into convolve_stamp callback

3. **Delta basis convolve_stamp already exists:**
   - Just needs to be called from fillStamp() through interface

4. **Future bases (PCA, Fourier):**
   - Implement convolve_stamp callback
   - No other changes needed to fillStamp()

## Conclusion

The architecture is **functionally sound but aesthetically incomplete**. The pluggable basis system works well for initialization and kernel evaluation dispatch, but the stamp convolution loop is a missed application of the pattern. 

This is acceptable because:
- ✓ It works correctly (verified in code review)
- ✓ It's well-documented and commented
- ✓ Adding delta basis didn't break anything
- ⚠ Extending with new bases requires fillStamp() modification
- ⚠ The unused convolve_stamp callback is confusing

Consider this refactoring as a nice-to-have for future basis implementations.
