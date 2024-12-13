struct SpectralEstimate{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	power::N
	function SpectralEstimate(freq, power)
        P = checkfreqdomaininputs(freq, power)
        new{P, typeof(freq), typeof(power)}(freq, power)
    end
end
getfreq(est::SpectralEstimate) = est.freq
getestimate(est::SpectralEstimate) = est.power

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
If data is a single process, and nfreq = (n_1,...,n_D), the output is a named tuple with
- freq: a tuple of array of frequencies in each dimension (note this is ordered in increasing frequency).
- power: a n_1 x ... x n_D array.

If the data has P processes, and nfreq = (n_1,...,n_D), the output is a named tuple with
- freq: a tuple of array of frequencies in each dimension (note this is ordered in increasing frequency).
- power: the P x P x n_1 x ... x n_D array.

"""
function multitaper_estimate(
	data::Union{GeoTable, PointSet},
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
	data::NTuple{P, Union{GeoTable, PointSet}},
	region;
	nfreq,
	fmax,
	tapers,
	mean_method = DefaultMean(),
) where {P}
	data, dim = check_spatial_data(data)
	mean_method = check_mean_method(mean_method, data)
	J_n = tapered_dft(data, tapers, nfreq, fmax, region, mean_method)
	freq = make_freq(nfreq, fmax, dim)
	power = dft2spectralmatrix(J_n)
	return SpectralEstimate(freq, power)
end

function preallocate_spectralmatrix(J_n::NTuple{1, Array{T, N}}) where {T, N}
	return Array{T, N - 1}(undef, size(J_n[1])[1:end-1])
end
function preallocate_spectralmatrix(J_n::NTuple{P, Array{T, N}}) where {P, T, N}
	return Array{SMatrix{P, P, T, P * P}, N - 1}(undef, size(J_n[1])[1:end-1])
end
function dft2spectralmatrix(J_n::NTuple{P, Array{T, N}}) where {P, T, N}
	S_mat = preallocate_spectralmatrix(J_n)
	dft2spectralmatrix!(S_mat, J_n)
	return S_mat
end

function dft2spectralmatrix!(
	S_mat::Array{T, D},
	J_n::NTuple{1, Array{T, N}},) where {T, N, D}
    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(abs2, @view J_n[1][i, :])
    end
    return S_mat
end

function dft2spectralmatrix!(
	S_mat::Array{SMatrix{P, P, T, L}, D},
	J_n::NTuple{P, Array{T, N}},
) where {P, T, N, L, D}
	# at this point J_n is a P-tuple of DFTs of dimension n_1 x ... x n_D x M
	# we want to return a n_1 x ... x n_D array of static P x P matrices (TODO maybe even hermitian symmetric ones later)
	@assert all(size(S_mat) == size(J)[1:end-1] for J in J_n) "S_mat should have the same size as the first N-1 dimensions of each J_n"
	for i in CartesianIndices(S_mat)
		S_mat[i] = mean(
			outerprod(SVector(ntuple(j -> J_n[j][i, m], Val{P}()))) for
			m in axes(J_n[1], N)
		)
	end
	return S_mat
end
outerprod(x) = x * x'

dft2spectralmatrix(J_n) = mapslices(x -> spectral_matrix(x), J_n, dims = (1, 2))

make_freq(nfreq::Int, fmax::Number, dim::Int) =
	ntuple(d -> choose_freq_1d(nfreq, fmax), dim)
function make_freq(nfreq, fmax, dim::Int)
	freq = choose_freq_1d.(nfreq, fmax)
	@assert length(freq) == dim "error in passing function, dim should be equal to length of freq"
	return freq
end
