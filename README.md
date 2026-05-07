# RegisterMismatchCommon.jl

[![CI](https://github.com/HolyLab/RegisterMismatchCommon.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/HolyLab/RegisterMismatchCommon.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/HolyLab/RegisterMismatchCommon.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HolyLab/RegisterMismatchCommon.jl)

RegisterMismatchCommon provides the shared types, utilities, and aperture-based
workflow helpers used for image-registration mismatch computation in the HolyLab
ecosystem. Concrete `mismatch` implementations live in downstream packages:

- [`RegisterMismatch`](https://github.com/HolyLab/RegisterMismatch.jl) — CPU/FFTW
- [`RegisterMismatchCuda`](https://github.com/HolyLab/RegisterMismatchCuda.jl) — GPU/CUFFT

In practice you will `using RegisterMismatch` (or its CUDA variant) rather than
this package directly. RegisterMismatchCommon is an explicit dependency for code
that needs to type-annotate or dispatch on `MismatchArray` / `NumDenom` without
caring which backend is loaded.

## Installation

This package lives in the [HolyLab registry](https://github.com/HolyLab/HolyLabRegistry).
Add the registry once, then install normally:

```julia
using Pkg
pkg"registry add General https://github.com/HolyLab/HolyLabRegistry.git"
Pkg.add("RegisterMismatchCommon")
```

## Key concepts

### `NumDenom` and `MismatchArray`

Mismatch is stored as a ratio: a *numerator* (sum of squared pixel differences)
over a *denominator* (normalization factor). Keeping them separate lets
downstream code combine apertures, apply thresholds, and interpolate without
losing information. `NumDenom{T}` is a two-field struct from
[`RegisterCore`](https://github.com/HolyLab/RegisterCore.jl), and a
`MismatchArray` is a `CenterIndexedArray` of `NumDenom` values indexed from
`-maxshift` to `+maxshift` along each dimension.

### Aperture workflow

Rather than computing a single mismatch over the whole image, the aperture
workflow divides the image into overlapping sub-regions (*apertures*) and
computes one `MismatchArray` per aperture. This yields a spatially-resolved
shift field, which is the starting point for non-rigid registration.

The typical sequence is:

```
aperture_grid         →   center coordinates for each aperture
allocate_mmarrays     →   pre-allocated output array
mismatch_apertures    →   fill the output  (needs a backend loaded)
correctbias!          →   remove pixel-bias artifacts
truncatenoise!        →   zero out low-signal entries
```

### `mismatch` requires a backend

`mismatch` and functions that call it (`mismatch_apertures`, `register_translate`)
are stubs defined here but not implemented. Loading `RegisterMismatch` or
`RegisterMismatchCuda` extends them with a concrete method. Calling them without
a backend will throw a `MethodError`.

## Usage examples

### Aperture grid

```julia
using RegisterMismatchCommon

# 64×64 image, 4×4 grid of apertures
ag = aperture_grid((64, 64), (4, 4))
size(ag)                      # (4, 4)
ag[1, 1]                      # (1.0, 1.0)   — top-left corner
ag[4, 4]                      # (64.0, 64.0) — bottom-right corner
```

### Aperture width

```julia
img = zeros(Float32, 64, 64)
default_aperture_width(img, (4, 4))   # (21.0, 21.0)
```

### Allocating output arrays

```julia
# Pre-allocate a 4×4 grid of MismatchArrays, each covering ±5 pixels
mms = allocate_mmarrays(Float32, (4, 4), (5, 5))
size(mms)        # (4, 4)
size(mms[1, 1])  # (11, 11)  — 2*5+1 per dimension
```

### Zero-shift mismatch (no backend needed)

`mismatch0` measures the mismatch at zero shift directly, without FFTs:

```julia
fixed  = [1.0 2.0; 3.0 4.0]
moving = [1.0 2.0; 3.0 4.0]
mismatch0(fixed, moving)          # NumDenom(0.0, 60.0) — perfect match

moving2 = [2.0 3.0; 4.0 5.0]
mm0 = mismatch0(fixed, moving2)   # NumDenom(4.0, 84.0)
mm0.num / mm0.denom               # ≈ 0.048 — normalized mismatch
```

### Padding mismatched arrays

```julia
a = [1.0 2.0; 3.0 4.0]           # 2×2
b = [5.0 6.0 7.0; 8.0 9.0 10.0]  # 2×3
ap, bp = nanpad(a, b)
# ap is 2×3, padded with NaN in the third column
# bp is returned unchanged
```

### Post-processing

```julia
# Zero out entries whose denominator is too small to be reliable
truncatenoise!(mms, 0.5f0)

# Impute zero-shift entries that are corrupted by pixel-bias
correctbias!(mms)
```

### Full-image translation (with a backend)

```julia
using RegisterMismatch   # or RegisterMismatchCuda

fixed  = rand(Float32, 64, 64)
moving = rand(Float32, 64, 64)

shift = register_translate(fixed, moving, (10, 10))
# CartesianIndex of best integer shift
```

## API reference

### Core computation

| Function | Description |
|---|---|
| `mismatch` | Full-array mismatch (requires backend) |
| `mismatch_apertures` | Aperture-wise mismatch on a grid (requires backend) |
| `mismatch0` | Zero-shift mismatch without FFTs |
| `register_translate` | Best integer translation (requires backend) |

### Aperture workflow

| Function | Description |
|---|---|
| `aperture_grid` | Uniformly-spaced grid of aperture center coordinates |
| `allocate_mmarrays` | Pre-allocate an array of `MismatchArray`s |
| `default_aperture_width` | Compute aperture width for a given grid |
| `aperture_range` | `UnitRange` indices for one aperture |
| `each_point` | Iterate over aperture centers in any layout |

### Post-processing

| Function | Description |
|---|---|
| `correctbias!` | Impute pixel-bias-corrupted zero-shift entries |
| `truncatenoise!` | Zero out low-denominator (unreliable) entries |

### Utilities

| Function | Description |
|---|---|
| `nanpad` | Pad two arrays to the same size with `NaN` |
| `padsize` | FFT-friendly padded size |
| `padranges` | Padded index ranges for FFT cross-correlation |
| `checksize_maxshift` | Validate a mismatch array's size |
| `assertsamesize` | Throw if two arrays differ in size |
| `tovec` | Convert a tuple to a `Vector` |
| `shiftrange` | Shift a range by a scalar |
| `set_FFTPROD` | Set the allowed FFT prime factors |

### Types

| Name | Description |
|---|---|
| `DimsLike` | `Union{AbstractVector{Int}, Dims}` — dimension-size argument |
| `WidthLike` | `Union{AbstractVector, Tuple}` — aperture-width argument |
