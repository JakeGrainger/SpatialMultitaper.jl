struct SpectralEstimate{D,F,P,N,T<:Union{Int, Nothing}} <: AnisotropicEstimate{D,P}
    freq::NTuple{D,F}
    power::N
    ntapers::T
    function SpectralEstimate(freq::NTuple{D,F}, power, ntapers) where {D,F}
        P = checkinputs(freq, power)
        new{D,F,P,typeof(power), typeof(ntapers)}(freq, power, ntapers)
    end
end
getargument(est::SpectralEstimate) = est.freq
getestimate(est::SpectralEstimate) = est.power
getextrafields(est::SpectralEstimate) = (est.ntapers,)

"""
	multitaper_estimate(data, region; nfreq, fmax, tapers, mean_method = DefaultMean())

Computes the multitaper spectral estimate from a tapered dft.

# Arguments:
- `data`: The data to estimate the spectrum from.
- `region`: The region to estimate the spectrum from.
- `nfreq::NTuple{D,Int}`: The number of frequencies in each dimension.
- `fmax::NTuple{D,Real}`: The maximum frequency in each dimension.
- `tapers`: A tuple of taper functions.
- `mean_method`: The method to estimate the mean. Default is `DefaultMean`, but can use `KnownMean(x)` to specify this if known.

# Output:
The output is a `SpectralEstimate` object, with a `freq` and `power` field.
Say `nfreq = (n_1,...,n_D)`, then `freq` is a `D` dimensional `NTuple` of arrays of frequencies in each dimension.
Power can be in the following forms:
- If `data` is a single process, then `power` is a `n_1 x ... x n_D` array.
- If `data` has `P` processes and is passed as an `NTuple`, then `power` is an `n_1 x ... x n_D` array of `SMatrix{P,P,T,P*P}`.
- If `data` is a vector of `P` processes, then `power` is a `P x P x n_1 x ... x n_D` array.

In any case, indexing into a `SpectralEstimate` will provide a `SpectralEstimate` object with the same `freq` and a subset of the `power` array corresponding to the processes in the chosen index.
"""
function multitaper_estimate(
    data::Union{GeoTable,PointSet},
    region;
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    return multitaper_estimate(
        (data,),
        region,
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
    )
end
function multitaper_estimate(
    data,
    region;
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    mask.(data, Ref(region))
    data, dim = check_spatial_data(data)
    mean_method = check_mean_method(mean_method, data)
    J_n = tapered_dft(data, tapers, nfreq, fmax, region, mean_method)
    freq = make_freq(nfreq, fmax, dim)
    power = dft2spectralmatrix(J_n)
    return SpectralEstimate(freq, power, length(tapers))
end

"""
    dft2spectralmatrix(J_n::AbstractArray)

Computes the spectral matrix from the DFTs, assuming that the DFTs are stored as one large array which is P x M x n_1 x ... x n_D.
"""
dft2spectralmatrix(J_n::AbstractArray) =
    mapslices(x -> spectral_matrix(x), J_n, dims = (1, 2))


"""
    dft2spectralmatrix(J_n::NTuple{P, AbstractArray{T, N}}) where {P, T, N}

Computes the spectral matrix from the DFTs, assuming that the DFTs are stored as a tuple of P arrays of size n_1 x ... x n_D x M.
"""
function dft2spectralmatrix(J_n::NTuple{P,AbstractArray{T,N}}) where {P,T,N}
    S_mat = preallocate_spectralmatrix(J_n)
    dft2spectralmatrix!(S_mat, J_n)
    return S_mat
end


function preallocate_spectralmatrix(J_n::NTuple{1,AbstractArray{T,N}}) where {T,N}
    return Array{T,N - 1}(undef, size(J_n[1])[1:end-1])
end
function preallocate_spectralmatrix(J_n::NTuple{P,AbstractArray{T,N}}) where {P,T,N}
    return Array{SMatrix{P,P,T,P * P},N - 1}(undef, size(J_n[1])[1:end-1])
end

function dft2spectralmatrix!(S_mat::Array{T,D}, J_n::NTuple{1,Array{T,N}}) where {T,N,D}
    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(abs2, @view J_n[1][i, :])
    end
    return S_mat
end

function dft2spectralmatrix!(
    S_mat::Array{SMatrix{P,P,T,L},D},
    J_n::NTuple{P,AbstractArray{T,N}},
) where {P,T,N,L,D}
    # at this point J_n is a P-tuple of DFTs of dimension n_1 x ... x n_D x M
    # we want to return a n_1 x ... x n_D array of static P x P matrices (TODO maybe even hermitian symmetric ones later)
    @assert all(size(S_mat) == size(J)[1:end-1] for J in J_n) "S_mat should have the same size as the first N-1 dimensions of each J_n"
    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(
            spectral_matrix(SVector(ntuple(j -> J_n[j][i, m], Val{P}()))) for
            m in axes(J_n[1], N)
        )
    end
    return S_mat
end

spectral_matrix(x::AbstractVector) = x * x'
spectral_matrix(x::AbstractMatrix) = (x * x') ./ size(x, 2)


make_freq(nfreq::Int, fmax::Number, dim::Int) =
    ntuple(d -> choose_freq_1d(nfreq, fmax), dim)
function make_freq(nfreq, fmax, dim::Int)
    freq = choose_freq_1d.(nfreq, fmax)
    @assert length(freq) == dim "error in passing function, dim should be equal to length of freq (which is a tuple of frequencies)"
    return freq
end
