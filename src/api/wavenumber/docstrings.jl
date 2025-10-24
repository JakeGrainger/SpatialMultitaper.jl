"""
    spectra(data; nk, kmax, dk, tapers, nw, mean_method) -> Spectra

Compute the multitaper spectral estimate from a tapered DFT.

# Arguments
- `data`: The data to estimate the spectrum from, of type `SpatialData`

# Keywords
- `nk = nothing`: The number of wavenumbers in each
    dimension. Can be scalar (applied uniformly) or tuple specifying each dimension.
- `kmax = nothing`: The maximum wavenumber in each dimension. Can be
    scalar (applied uniformly) or tuple specifying each dimension.
- `dk = nothing`: The wavenumber spacing in each
    dimension. Can be scalar (applied uniformly) or tuple specifying each dimension.
- `tapers = nothing`: A tuple of taper functions.
    If `nothing`, tapers will be generated using `nw`.
- `nw = 3`: The space-bandwidth product for generating tapers when `tapers` is not
    provided. Represents side length times half-bandwidth (`nw = L * W`).
- `mean_method::MeanEstimationMethod = DefaultMean()`: The method to estimate the mean.

You only need to specify two of the three parameters `nk`, `kmax`, and `dk`. If one of
`nk` and `kmax` is specified, `dk` will be set to a default based on the region. Parameters
can be mixed as scalars and tuples/vectors, e.g. `nk=100` and `kmax=(0.5, 1.0)` for 2D data.
`nk` values must be positive integers (Real values are rounded up), `kmax` and `dk` must be
positive Real numbers.

# Returns
- `Spectra`: A spectral estimate object with `wavenumber` and `power` fields:
  - `wavenumber`: D-dimensional `NTuple` of wavenumber arrays for each dimension
  - `power`: Power spectral density with shape depending on data type:
    - Single process: `n_1 × ... × n_D` array
    - `NTuple{P}` data: `n_1 × ... × n_D` array of `SMatrix{P, P}`
    - Vector of P processes: `P × P × n_1 × ... × n_D` array

# Notes
- Indexing into a `Spectra` object with two indices indexes into the process dimensions,
    e.g. `S[1, 2]` gives the cross-spectrum between processes 1 and 2.
- Use `KnownMean(x)` to specify a known mean value
- You can also pass `data` and `region` as arguments which will first be passed to
    `spatial_data` to construct a `SpatialData` object.
- For irregular regions, `nw` represents the radius of a circle in wavenumber space on
    which tapers are concentrated. For more fine control, construct your own tapers.

# Examples
```julia
# Basic spectral estimation with automatic taper generation
spec = spectra(data; kmax = 0.5, nw = 3)

# Using custom wavenumber parameters
spec = spectra(data; nk = 64, kmax = (0.5, 1.0), nw = 4)

# With known mean
spec = spectra(data; kmax = 0.3, mean_method = KnownMean(0.0))
```
"""
spectra

"""
    partial_spectra(spectrum::Spectra{MarginalTrait}) -> Spectra{PartialTrait}
    partial_spectra(spectrum::RotationalSpectra{MarginalTrait}) -> RotationalEstimate{PartialTrait}
    partial_spectra(data, region; kwargs...) -> Spectra{PartialTrait}
    partial_spectra(data::SpatialData; kwargs...) -> Spectra{PartialTrait}

Compute partial spectral estimates from marginal spectral estimates or directly from data.

Partial spectra remove the linear influence of all other processes, providing a measure of
the direct relationship between process pairs. Unlike marginal spectra which show total
power including indirect effects, partial spectra reveal only the direct linear
relationships after removing the influence of all other processes. The computation involves
matrix inversion of the spectral matrix and finite-sample bias correction based on the
number of tapers used in the original estimation.

# Arguments
- `spectrum::Spectra{MarginalTrait}`: A marginal spectral estimate to convert
- `spectrum::RotationalSpectra{MarginalTrait}`: A rotational marginal spectral estimate
- `data`: Spatial data for direct partial spectra computation
- `region::Meshes.Geometry`: Spatial region for direct computation

# Keywords
When computing directly from data, all keywords from [`spectra`](@ref) are supported:
- `nk`: Number of wavenumbers in each dimension
- `kmax`: Maximum wavenumber in each dimension
- `dk`: Wavenumber spacing in each dimension
- `tapers`: Taper functions to use
- `nw = 3`: Space-bandwidth product for taper generation
- `mean_method = DefaultMean()`: Method for mean estimation

# Returns
- `Spectra{PartialTrait}`: Partial spectral estimate with same wavenumber grid as input
- `RotationalEstimate{PartialTrait}`: For rotational input spectra

The returned partial spectra have the same spatial dimensions and wavenumber grid as the
input, but represent direct relationships between processes rather than total power.

# Throws
- `ArgumentError`: If the spectrum does not have equal input and output process sets
    (i.e., the spectral matrix is not square). Partial spectra require square spectral
    matrices, so you shouldn't have subsetted before calling partial spectra.

# Mathematical Details
For a spectral matrix `f`, the partial spectrum `P` is computed through matrix inversion:

**Basic Formula:**
- Diagonal elements: `Pᵢᵢ = 1/Gᵢᵢ`
- Off-diagonal elements: `Pᵢⱼ = -Gᵢⱼ/(GᵢᵢGⱼⱼ - |Gᵢⱼ|²)`

where `G = f⁻¹` is the inverse spectral matrix.

**Bias Correction:**
Finite-sample bias correction is automatically applied using the number of tapers:
- Correction factor: `M/(M - Q + xᵢⱼ)` where:
  -`M = number of tapers from the original spectral estimate
  - `Q` = number of processes
  - `xᵢⱼ` = 1 if i=j (diagonal), 2 if i≠j (off-diagonal)

**Special Cases:**
- Single process: Returns the original spectral value unchanged
- Rotational spectra: Bias correction is handled differently (not fully implemented)

# Notes
- Partial spectra are the spectral domain equivalent of partial correlation
- Values can be complex and may have larger magnitudes than marginal spectra
- Diagonal elements represent the "partial power" of each process
- Off-diagonal elements show direct cross-relationships
- Use [`partial_spectra_uncorrected`](@ref) to skip bias correction
- For rotational estimates, bias correction is not currently fully implemented

# Examples
```julia
# Compute partial spectra from existing marginal estimates
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
partial_spec = partial_spectra(marginal_spec)

# Direct computation from data and region
partial_spec = partial_spectra(data, region; nk = (32, 32), kmax = (0.5, 0.5), nw = 4)

# Direct computation from SpatialData object
spatial_data_obj = spatial_data(data, region)
partial_spec = partial_spectra(spatial_data_obj; kmax = 0.3, tapers = my_tapers)

# Access partial relationships between processes
direct_coupling = partial_spec[1, 2]  # Direct coupling between processes 1 and 2
partial_power = partial_spec[1, 1]    # Partial power of process 1
```

See also: [`spectra`](@ref), [`partial_coherence`](@ref), [`partial_spectra_uncorrected`](@ref)
"""
partial_spectra

"""
    coherence(spectrum::NormalOrRotationalSpectra{E}) where {E} -> Coherence{E}
    coherence(data, region; kwargs...) -> Coherence
    coherence(data::SpatialData; kwargs...) -> Coherence

Compute coherence from a spectral estimate or directly from spatial data.

The coherence function computes the normalized cross-spectral density, providing a measure
of linear dependence between processes as a function of wavenumber. For a cross-spectral
matrix S, the coherence γᵢⱼ is computed as γᵢⱼ = Sᵢⱼ / √(Sᵢᵢ * Sⱼⱼ).

# Arguments
- `spectrum::NormalOrRotationalSpectra{E}`: A `Spectra` or `RotationalSpectra` estimate
    with multiple processes
- `data`: Spatial data for direct coherence computation
- `region::Meshes.Geometry`: Spatial region for direct computation

# Keywords
When computing directly from data, all keywords from [`spectra`](@ref) are supported:
- `nk`: Number of wavenumbers in each dimension
- `kmax`: Maximum wavenumber in each dimension
- `dk`: Wavenumber spacing in each dimension
- `tapers`: Taper functions to use
- `nw = 3`: Space-bandwidth product for taper generation
- `mean_method = DefaultMean()`: Method for mean estimation

# Returns
- `Coherence{E}`: A coherence estimate object containing:
  - `wavenumber`: Wavenumber grid matching the input spectrum
  - `coherence`: Coherence estimates with the same spatial dimensions as input
  - `processinformation`: Information about the analyzed processes
  - `estimationinformation`: Details about the estimation procedure

# Throws
- `ArgumentError`: If the spectrum does not have equal input and output process sets
    (i.e., the spectral matrix is not square)

# Notes
- Coherence values are complex numbers with magnitude ≤ 1
- Diagonal elements are always 1 (perfect self-coherence)
- For single-process data, returns scalar coherence of 1
- Use [`magnitude_coherence`](@ref) or [`magnitude_squared_coherence`](@ref) for
    magnitude-only results

# Examples
```julia
# Compute coherence from existing spectral estimate
spec = spectra(data; kmax = 0.5, nw = 3)
coh = coherence(spec)

# Direct computation from data and region
coh = coherence(data, region; nk = (32, 32), kmax = (0.5, 0.5), nw = 4)

# Direct computation from SpatialData object
spatial_data_obj = spatial_data(data, region)
coh = coherence(spatial_data_obj; kmax = 0.3, tapers = my_tapers)

# Access coherence between processes 1 and 2
cross_coherence = coh[1, 2]
```
"""
coherence

"""
    partial_coherence(spectrum) -> Coherence{PartialTrait}
    partial_coherence(data, region; kwargs...) -> Coherence{PartialTrait}
    partial_coherence(data::SpatialData; kwargs...) -> Coherence{PartialTrait}

Compute partial coherence from a spectral estimate or directly from spatial data.

See [`coherence`](@ref) for details on arguments, keywords, and return types.
"""
partial_coherence
