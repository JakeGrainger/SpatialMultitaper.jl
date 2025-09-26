struct Spectra{E, D, P, Q, A, T, IP, IE} <: AnisotropicEstimate{E, D, P, Q}
    freq::NTuple{D, A}
    power::T
    processinformation::IP
    estimationinformation::IE
    function Spectra{E}(freq::NTuple{D, A}, power::T, processinfo::IP,
            estimationinfo::IE) where {E <: EstimateTrait, D, A, T, IP, IE}
        P, Q = checkinputs(freq, power, processinfo)
        new{E, D, P, Q, A, T, IP, IE}(freq, power, processinfo, estimationinfo)
    end
end

const RotationalSpectra{E, D, P, Q, S <: Spectra} = RotationalEstimate{E, D, P, Q, S}
const NormalOrRotationalSpectra{E} = Union{Spectra{E}, RotationalSpectra{E}}
getargument(est::Spectra) = est.freq
getestimate(est::Spectra) = est.power

"""
	spectra(data, region; nfreq, fmax, tapers, mean_method = DefaultMean())

Computes the multitaper spectral estimate from a tapered dft.

# Arguments:
- `data`: The data to estimate the spectrum from.
- `region`: The region to estimate the spectrum from.
- `nfreq::NTuple{D,Int}`: The number of frequencies in each dimension.
- `fmax::NTuple{D,Real}`: The maximum frequency in each dimension.
- `tapers`: A tuple of taper functions.
- `mean_method`: The method to estimate the mean. Default is `DefaultMean`.

# Output:
The output is a `Spectra` object, with a `freq` and `power` field.
Say `nfreq = (n_1,...,n_D)`, then `freq` is a `D` dimensional `NTuple` of arrays of
frequencies in each dimension.
Power can be in the following forms:
- If `data` is a single process, then `power` is a `n_1 x ... x n_D` array.
- If `data` is an `NTuple{P}`, then `power` is a `n_1 x ... x n_D` array of `SMatrix{P, P}`.
- If `data` is a vector of `P` processes, then `power` is a `P x P x n_1 x ... x n_D` array.

In any case, indexing into a `Spectra` will provide a `Spectra` object with the same `freq`
and a subset of the `power` array corresponding to the processes in the chosen index.

`KnownMean(x)` can be used to specify the mean when it is known.
"""
function spectra(data, region::Meshes.Geometry; kwargs...)
    return spectra(spatial_data(data, region); kwargs...)
end
function spectra(data::SpatialData; nfreq, fmax, tapers,
        mean_method::MeanEstimationMethod = DefaultMean())
    freq = make_freq(nfreq, fmax, embeddim(data))
    J_n = tapered_dft(data, tapers, nfreq, fmax, mean_method)
    power = dft2spectralmatrix(J_n, Val{ncol(data)}())

    processinformation = ProcessInformation(data; mean_method = mean_method)
    estimationinformation = EstimationInformation(length(tapers))

    return Spectra{MarginalTrait}(
        freq, power, processinformation, estimationinformation)
end

multitaper_estimate = spectra

"""
    dft2spectralmatrix(J_n::AbstractArray, Val{P})

Computes the spectral matrix from the DFTs, assuming that the DFTs are stored as one large
array which is P x M x n_1 x ... x n_D.
"""
function dft2spectralmatrix(J_n::AbstractArray, ::Val{P}) where {P}
    @argcheck size(J_n, 1) == P
    mapslices(x -> spectral_matrix(x), J_n, dims = (1, 2))
end

"""
    dft2spectralmatrix(J_n::AbstractArray, Val{1})

Computes the spectral matrix from the DFTs, assuming that we have a single DFT stored as one
large which is n_1 x ... x n_D x M.
"""
function dft2spectralmatrix(J_n::AbstractArray, ::Val{1})
    reshape(mapslices(x -> mean(abs2, x), J_n, dims = ndims(J_n)), size(J_n)[1:(end - 1)])
end

"""
    dft2spectralmatrix(J_n::NTuple{P, AbstractArray{T, N}}, Val{P}) where {P, T, N}

Computes the spectral matrix from the DFTs, assuming that the DFTs are stored as a tuple of
P arrays of size n_1 x ... x n_D x M.
"""
function dft2spectralmatrix(J_n::NTuple{P, AbstractArray{T, N}}, ::Val{P}) where {P, T, N}
    S_mat = preallocate_spectralmatrix(J_n)
    dft2spectralmatrix!(S_mat, J_n)
    return S_mat
end

function preallocate_spectralmatrix(J_n::NTuple{P, AbstractArray{T, N}}) where {P, T, N}
    return Array{SMatrix{P, P, T, P * P}, N - 1}(undef, size(J_n[1])[1:(end - 1)])
end

function dft2spectralmatrix!(
        S_mat::Array{T, D}, J_n::NTuple{1, Array{T, N}}) where {T, N, D}
    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(abs2, @view J_n[1][i, :])
    end
    return S_mat
end

function dft2spectralmatrix!(
        S_mat::Array{<:SMatrix{P, P, T}},
        J_n::NTuple{P, AbstractArray{T, N}}
) where {P, T, N}
    # at this point J_n is a P-tuple of DFTs of dimension n_1 x ... x n_D x M
    # we want to return a n_1 x ... x n_D array of static P x P matrices (TODO maybe even hermitian symmetric ones later)
    @assert all(size(S_mat) == size(J)[1:(end - 1)] for J in J_n) "S_mat should have the same size as the first N-1 dimensions of each J_n"
    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(spectral_matrix(SVector(ntuple(j -> J_n[j][i, m], Val{P}())))
        for m in axes(J_n[1], N))
    end
    return S_mat
end

spectral_matrix(x::AbstractVector) = x * x'
spectral_matrix(x::AbstractMatrix) = (x * x') ./ size(x, 2)

function make_freq(nfreq::Int, fmax::Number, dim::Int)
    ntuple(d -> choose_freq_1d(nfreq, fmax), dim)
end
function make_freq(nfreq, fmax, dim::Int)
    freq = choose_freq_1d.(nfreq, fmax)
    @assert length(freq)==dim "error in passing function, dim should be equal to length of freq (which is a tuple of frequencies)"
    return freq
end
