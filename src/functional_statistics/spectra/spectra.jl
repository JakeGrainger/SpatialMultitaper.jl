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

struct Spectra{E, D, A, T, IP, IE} <: AnisotropicEstimate{E, D}
    wavenumber::NTuple{D, A}
    power::T
    processinformation::IP
    estimationinformation::IE
    function Spectra{E}(wavenumber::NTuple{D, A}, power::T, processinfo::IP,
            estimationinfo::IE) where {E <: EstimateTrait, D, A, T, IP, IE}
        checkinputs(wavenumber, power, processinfo)
        return new{E, D, A, T, IP, IE}(wavenumber, power, processinfo, estimationinfo)
    end
end

const RotationalSpectra{E, D, S <: Spectra} = RotationalEstimate{E, D, S}
const NormalOrRotationalSpectra{E} = Union{Spectra{E}, RotationalSpectra{E}}
getshortestimatename(::Type{<:Spectra}) = "f"
get_evaluation_points(est::Spectra) = est.wavenumber
get_estimates(est::Spectra) = est.power

function spectra(data, region::Meshes.Geometry; kwargs...)::Spectra
    return spectra(spatial_data(data, region); kwargs...)
end

function spectra(data::SpatialData; nk = nothing, kmax = nothing, dk = nothing,
        tapers = nothing, nw = 3, mean_method::MeanEstimationMethod = DefaultMean())::Spectra
    tapers = _validate_tapers(tapers, getregion(data), nw)
    all_mem = preallocate_spectra(data; nk = nk, kmax = kmax, dk = dk, tapers = tapers)
    return spectra!(all_mem, data; nk = nk, kmax = kmax, dk = dk,
        tapers = tapers, mean_method = mean_method)
end

function preallocate_spectra(
        data::SpatialData; nk = nothing, kmax, dk = nothing, tapers = nothing, nw = 3)
    _nk, _kmax = _validate_wavenumber_params(nk, kmax, dk, data)
    tapers = _validate_tapers(tapers, getregion(data), nw)
    mem = preallocate_tapered_dft(data, tapers, _nk, _kmax)
    power = _preallocate_spectral_matrix(data, _nk)
    return (mem, power)
end

function spectra!(
        all_mem, data::SpatialData; nk = nothing, kmax, dk = nothing,
        tapers = nothing, nw = 3, mean_method::MeanEstimationMethod = DefaultMean())::Spectra
    _nk, _kmax = _validate_wavenumber_params(nk, kmax, dk, data)
    wavenumber = _make_wavenumber_grid(_nk, _kmax)
    tapers = _validate_tapers(tapers, getregion(data), nw)

    mem = all_mem[1]
    power = all_mem[2]

    J_n = tapered_dft!(mem, data, tapers, _nk, _kmax, mean_method)
    _dft_to_spectral_matrix!(power, J_n, process_trait(data))

    process_info = ProcessInformation(data; mean_method = mean_method)
    estimation_info = EstimationInformation(length(tapers))

    return Spectra{MarginalTrait}(wavenumber, power, process_info, estimation_info)
end

# Alias for backwards compatibility
const multitaper_estimate = spectra

"""
    _dft_to_spectral_matrix(data, J_n, nk)

Compute the spectral matrix from DFTs for different process types.

# Arguments
- `J_n`: The tapered DFTs, with shape depending on the process type:
    - For `SingleProcessTrait`: `n₁ × ... × n_D × M` array
    - For `MultipleVectorTrait`: `P × M × n₁ × ... × n_D` array
    - For `MultipleTupleTrait`: A `n₁ × ... × n_D × M` array of `SVector`s of length `P`
- Trait type: Specifies the process structure.

# Returns
- The spectral matrix or array of spectral matrices, with shape depending on the input.

# Notes
- For single process: returns an array of power spectral densities.
- For multiple processes (tuple): returns an array of spectral matrices.
- For multiple vector processes: returns a D+2 array where each slice is a spectral matrix.
"""
function _dft_to_spectral_matrix(data, J_n, nk)
    S_mat = _preallocate_spectral_matrix(data, nk)
    _dft_to_spectral_matrix!(S_mat, J_n, process_trait(data))
    return S_mat
end

function _preallocate_spectral_matrix(data::MultipleSpatialDataVec, nk)
    return zeros(ComplexF64, (ncol(data), ncol(data), nk...)) # TODO: generalize type
end

function _dft_to_spectral_matrix!(
        S_mat::AbstractArray, J_n::AbstractArray, ::MultipleVectorTrait)
    for i in CartesianIndices(size(J_n)[3:end])
        @views _compute_spectral_matrix!(S_mat[:, :, i], J_n[:, :, i])
    end
    return S_mat
end

function _preallocate_spectral_matrix(::SingleProcessData, nk)
    return zeros(Float64, nk) # TODO: generalize type
end

function _dft_to_spectral_matrix!(
        S_mat::AbstractArray, J_n::AbstractArray, ::SingleProcessTrait)
    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(abs2, @view J_n[i, :])
    end
    return S_mat
end

function _preallocate_spectral_matrix(data::MultipleSpatialDataTuple{P}, nk) where {P}
    return zeros(SMatrix{P, P, ComplexF64, P * P}, nk) # TODO: generalize type
end

"""
    _dft_to_spectral_matrix!(
        S_mat::Array{<:SMatrix{P, P}}, J_n::AbstractArray{<:SVector{P}, N},
        ::MultipleTupleTrait) where {P, N}

Fill spectral matrix for multiple process case.

Each array in J_n has size n_1 × ... × n_D × M.
"""
function _dft_to_spectral_matrix!(
        S_mat::Array{<:SMatrix{P, P}}, J_n::AbstractArray{<:SVector{P}, N},
        ::MultipleTupleTrait) where {P, N}
    # Validate dimensions
    @argcheck size(S_mat) == size(J_n)[1:(end - 1)]

    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(_compute_spectral_matrix(J_n[i, m]) for m in axes(J_n, N))
    end
    return S_mat
end

"""
    _compute_spectral_matrix(x::AbstractVector)

Compute the outer product spectral matrix for a vector.
"""
_compute_spectral_matrix(x::AbstractVector) = x * x'

"""
    _compute_spectral_matrix(x::AbstractMatrix)

Compute the averaged spectral matrix for a matrix.
"""
_compute_spectral_matrix(x::SMatrix) = (x * x') ./ size(x, 2)

function _compute_spectral_matrix!(y, x::AbstractMatrix)
    for i in axes(y, 1), j in axes(y, 2)
        y[i, j] = @views mean(x[i, :] .* conj.(x[j, :]))
    end
end

"""
    _make_wavenumber_grid(nk, kmax)

Create a wavenumber grid with dimension-specific parameters.
"""
function _make_wavenumber_grid(nk, kmax)
    wavenumber = _choose_wavenumbers_1d.(nk, kmax)
    return wavenumber
end
