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
getargument(est::Spectra) = est.wavenumber
getestimate(est::Spectra) = est.power

"""
    spectra(data, region; nk, kmax, tapers, mean_method = DefaultMean())

Compute the multitaper spectral estimate from a tapered DFT.

# Arguments
- `data`: The data to estimate the spectrum from
- `region`: The region to estimate the spectrum from
- `nk`: The number of wavenumbers in each dimension
- `kmax`: The maximum wavenumber in each dimension
- `dk`: The wavenumber spacing in each dimension
- `tapers`: A tuple of taper functions
- `mean_method::MeanEstimationMethod`: The method to estimate the mean (default: `DefaultMean()`)

If one of `nk` and `kmax` is specified, `dk` will be set to a default based on the region.
Otherwise, you only need to specify two of the three parameters `nk`, `kmax`, and `dk`. They
can either be scalars (applied uniformly across all dimensions) or tuples specifying each
dimension. You can mix the two styles, e.g. `nk=100` and `kmax=(0.5, 1.0)` for 2D data.
`nk` must be an `Int` or tuple of `Int`s and should be positive, `kmax` and `dk` must be
`Real` or tuples of `Real`s and also positive.

# Returns
A `Spectra` object with `wavenumber` and `power` fields:
- `wavenumber`: D-dimensional `NTuple` of wavenumber arrays for each dimension
- `power`: Power spectral density in one of the following forms:
  - Single process: `n_1 × ... × n_D` array
  - `NTuple{P}` data: `n_1 × ... × n_D` array of `SMatrix{P, P}`
  - Vector of P processes: `P × P × n_1 × ... × n_D` array

# Notes
- Indexing into a `Spectra` object returns a subset with the same `wavenumber`
- Use `KnownMean(x)` to specify a known mean value

# Examples
```julia
spec = spectra(data, region, kmax=(0.5, 0.5), tapers=tapers)
```
"""
function spectra(data, region::Meshes.Geometry; kwargs...)::Spectra
    return spectra(spatial_data(data, region); kwargs...)
end

function spectra(data::SpatialData; nk = nothing, kmax, dk = default_dk(data, nk, kmax),
        tapers, mean_method::MeanEstimationMethod = DefaultMean())::Spectra
    _nk, _kmax = _validate_wavenumber_params(nk, kmax, dk, data)
    wavenumber = _make_wavenumber_grid(_nk, _kmax)
    J_n = tapered_dft(data, tapers, _nk, _kmax, mean_method)
    power = _dft_to_spectral_matrix(J_n, process_trait(data))

    process_info = ProcessInformation(data; mean_method = mean_method)
    estimation_info = EstimationInformation(length(tapers))

    return Spectra{MarginalTrait}(wavenumber, power, process_info, estimation_info)
end

# Alias for backwards compatibility
const multitaper_estimate = spectra

"""
    _dft_to_spectral_matrix(J_n::AbstractArray, ::SingleProcessTrait)
    _dft_to_spectral_matrix(J_n::AbstractArray, ::MultipleVectorTrait)
    _dft_to_spectral_matrix(J_n::NTuple{P, AbstractArray}, ::MultipleTupleTrait) where {P}

Compute the spectral matrix from DFTs for different process types.

# Arguments
- `J_n`: The tapered DFTs, with shape depending on the process type:
    - For `SingleProcessTrait`: `n₁ × ... × n_D × M` array (single process)
    - For `MultipleVectorTrait`: `P × M × n₁ × ... × n_D` array (multiple vector processes)
    - For `MultipleTupleTrait`: Tuple of `P` arrays, each of size `n₁ × ... × n_D × M` (multiple processes as tuple)
- Trait type: Specifies the process structure.

# Returns
- The spectral matrix or array of spectral matrices, with shape depending on the input.

# Notes
- For single process: returns an array of power spectral densities.
- For multiple processes (tuple): returns an array of spectral matrices.
- For multiple vector processes: returns a D+2 array where each slice is a spectral matrix.
"""
function _dft_to_spectral_matrix(J_n, trait)
    S_mat = _preallocate_spectral_matrix(J_n, trait)
    _fill_spectral_matrix!(S_mat, J_n, trait)
    return S_mat
end

function _preallocate_spectral_matrix(J_n::AbstractArray, ::MultipleVectorTrait)
    return zeros(eltype(J_n), (size(J_n, 1), size(J_n, 1), size(J_n)[3:end]...))
end

function _fill_spectral_matrix!(
        S_mat::AbstractArray, J_n::AbstractArray, ::MultipleVectorTrait)
    for i in CartesianIndices(size(J_n)[3:end])
        S_mat[:, :, i] = @views _compute_spectral_matrix(J_n[:, :, i])
    end
    return S_mat
end

function _preallocate_spectral_matrix(J_n::AbstractArray, ::SingleProcessTrait)
    return zeros(real(eltype(J_n)), size(J_n)[1:(end - 1)])
end

function _fill_spectral_matrix!(
        S_mat::AbstractArray, J_n::AbstractArray, ::SingleProcessTrait)
    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(abs2, @view J_n[i, :])
    end
    return S_mat
end

"""
    _preallocate_spectral_matrix(J_n::NTuple{P, AbstractArray{T, N}}) where {P, T, N}

Preallocate storage for the spectral matrix computation.
"""
function _preallocate_spectral_matrix(
        J_n::NTuple{P, AbstractArray{T, N}}, ::MultipleTupleTrait) where {P, T, N}
    output_dims = size(J_n[1])[1:(end - 1)]
    return Array{SMatrix{P, P, T, P * P}, N - 1}(undef, output_dims)
end

"""
    _fill_spectral_matrix!(S_mat::Array{T, D}, J_n::NTuple{1, Array{T, N}}) where {T, N, D}

Fill spectral matrix for single process case (P=1).
"""
function _fill_spectral_matrix!(S_mat::Array{T, D}, J_n::NTuple{1, Array{T, N}},
        ::MultipleTupleTrait) where {T, N, D}
    if !all(size(S_mat) == size(J)[1:(end - 1)] for J in J_n)
        throw(DimensionMismatch("S_mat dimensions must match first N-1 dimensions of each J_n array"))
    end
    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(abs2, @view J_n[1][i, :])
    end
    return S_mat
end

"""
    _fill_spectral_matrix!(S_mat::Array{<:SMatrix{P, P, T}}, J_n::NTuple{P, AbstractArray{T, N}}) where {P, T, N}

Fill spectral matrix for multiple process case.

Each array in J_n has size n_1 × ... × n_D × M.
"""
function _fill_spectral_matrix!(
        S_mat::Array{<:SMatrix{P, P, T}}, J_n::NTuple{P, AbstractArray{T, N}},
        ::MultipleTupleTrait) where {P, T, N}
    # Validate dimensions
    if !all(size(S_mat) == size(J)[1:(end - 1)] for J in J_n)
        throw(DimensionMismatch("S_mat dimensions must match first N-1 dimensions of each J_n array"))
    end

    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(_compute_spectral_matrix(SVector(ntuple(
                            j -> J_n[j][i, m], Val{P}())))
        for m in axes(J_n[1], N))
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
_compute_spectral_matrix(x::AbstractMatrix) = (x * x') ./ size(x, 2)

"""
    _make_wavenumber_grid(nk, kmax)

Create a wavenumber grid with dimension-specific parameters.
"""
function _make_wavenumber_grid(nk, kmax)
    wavenumber = _choose_wavenumbers_1d.(nk, kmax)
    return wavenumber
end
