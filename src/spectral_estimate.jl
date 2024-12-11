function check_spatial_data(data::NTuple{N, Union{GeoTable, PointSet}}) where {N}
	_getdims(x) = x isa PointSet ? embeddim(x) : embeddim(domain(x))
	dim = _getdims(first(data))
	@assert all(_getdims.(data) .== dim) "data should all be the same spatial dimension"
	copied_data = deepcopy(data)
	for proc in copied_data
		if proc isa GeoTable
			if length(values(proc)) > 1
				@warn "more than one random field provided to a geotable, currently we only process the first of these!"
			end
			if any(x -> abs(x) == Inf, values(proc)[1])
				error("Some fields have infinite values!")
			end
			if any(isnan, values(proc)[1])
				replace!(values(proc)[1], NaN => zero(eltype(values(proc)[1])))
			end
		end
	end
	return copied_data, dim
end

function check_mean_method(
	mean_method::MeanEstimationMethod,
	data::NTuple{N, Union{GeoTable, PointSet}},
) where {N}
	return ntuple(i -> mean_method, Val{N}())
end

function check_mean_method(
	mean_method::NTuple{P, MeanEstimationMethod},
	data::NTuple{N, Union{GeoTable, PointSet}},
) where {P, N}
	P === N ||
		throw(ArgumentError("Number of mean methods should match number of processes"))
	return mean_method
end

struct SpectralEstimate{P, F, N, J <: Union{Nothing, Vector{N}}}
	freq::F
	power::N
	power_jackknifed::J
	function SpectralEstimate(freq::F, power::N, power_jackknifed::J, ::Val{P}) where {P, F, N, J}
		
		new{P, F, N, J}(freq, power, power_jackknifed)
	end
end

"""
	multitaper_estimate(data, tapers; nfreq, fmax, region, jackknife=false, mean_method = DefaultMean())

Computes the multitaper spectral estimate from a tapered dft.

# Arguments:
- `data`: The data to estimate the spectrum from.
- `tapers`: A tuple of taper functions.
- `nfreq::NTuple{D,Int}`: The number of frequencies in each dimension.
- `fmax::NTuple{D,Real}`: The maximum frequency in each dimension.
- `region`: The region to estimate the spectrum from.
- `jackknife`: If true, jackknife (leave one out) replicates returned.
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
	jackknife = false,
	mean_method::MeanEstimationMethod = DefaultMean(),
)
	mt_est = multitaper_estimate(
		(data,),
		region,
		nfreq = nfreq,
		fmax = fmax,
		tapers = tapers,
		jackknife = jackknife,
		mean_method = mean_method,
	)
	if jackknife
		return SpectralEstimate(
			mt_est.freq,
			reshape(mt_est.power, size(mt_est.power)[3:end]),
			reshape.(mt_est.power_jackknifed, size(mt_est.power)[3:end]),
			Val{1}(),
		)
	else
		return SpectralEstimate(
			mt_est.freq,
			reshape(mt_est.power, size(mt_est.power)[3:end]),
			nothing,
			Val{1}(),
		)
	end
end
function multitaper_estimate(
	data::NTuple{P, Union{GeoTable, PointSet}},
	region;
	nfreq,
	fmax,
	tapers,
	jackknife = false,
	mean_method = DefaultMean(),
) where {P}
	data, dim = check_spatial_data(data)
	mean_method = check_mean_method(mean_method, data)
	J_n = tapered_dft(data, tapers, nfreq, fmax, region, mean_method)
	freq = make_freq(nfreq, fmax, dim)
    power = mapslices(x -> spectral_matrix(x), J_n, dims = (1, 2))
	if jackknife
		jk_weights = [make_jk_weight(size(J_n, 1), m) for m in axes(J_n, 1)] # weight vector which zeros out mth entry
		power_jackknifed = [
			mapslices(x -> spectral_matrix(x, jk_weights[m]), J_n, dims = (1, 2)) for
			m in axes(J_n, 1)
		]
		return SpectralEstimate(freq, power, power_jackknifed, Val{P}())
	else
		return SpectralEstimate(freq, power, nothing, Val{P}())
	end
end

function preallocate_spectralmatrix(J_n::NTuple{P, Array{T, N}}) where {P, T, N}
    return Array{SMatrix{P,P,T,P*P}, N-1}(undef, size(J_n[1])[1:end-1])
end
function dft2spectralmatrix(J_n::NTuple{P, Array{T, N}}) where {P, T, N}
    S_mat = preallocate_spectralmatrix(J_n)
    dft2spectralmatrix!(S_mat, J_n)
    return S_mat
end

function dft2spectralmatrix!(S_mat::Array{SMatrix{P,P,T,L}, D}, J_n::NTuple{P, Array{T, N}}) where {P, T, N, L, D}
    # at this point J_n is a P-tuple of DFTs of dimension n_1 x ... x n_D x M
    # we want to return a n_1 x ... x n_D array of static P x P matrices (maybe even hermitian symmetric ones later)
    @assert all(size(S_mat) == size(J)[1:end-1] for J in J_n) "S_mat should have the same size as the first N-1 dimensions of each J_n"
    for J in J_n
        for (i,sl) in enumerate(eachslice(J, dims = ntuple(i -> i, Val{N-1}())))
            J = SVector{P,T}(sl)
            S[i] = J*J'
        end
    end
    S_mat ./ P
    return S
end

dft2spectralmatrix(J_n) = mapslices(x -> spectral_matrix(x), J_n, dims = (1, 2))

function make_jk_weight(M, m)
	x = fill(1 / sqrt(M - 1), M)
	x[m] = 0
	return x
end

function spectral_matrix(x, weight = 1 / sqrt(size(x, 1)))
	xw = conj.(x) .* weight
	return xw' * xw
end

make_freq(nfreq::Int, fmax::Number, dim::Int) =
	ntuple(d -> choose_freq_1d(nfreq, fmax), dim)
function make_freq(nfreq, fmax, dim::Int)
	freq = choose_freq_1d.(nfreq, fmax)
	@assert length(freq) == dim "error in passing function, dim should be equal to length of freq"
	return freq
end
