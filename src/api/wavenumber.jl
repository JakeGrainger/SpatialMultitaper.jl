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

function spectra(data, region::Meshes.Geometry; kwargs...)::Spectra
    return spectra(spatial_data(data, region); kwargs...)
end

function spectra(data::SpatialData; kwargs...)::Spectra
    return compute(Spectra{MarginalTrait}, data; kwargs...)
end
