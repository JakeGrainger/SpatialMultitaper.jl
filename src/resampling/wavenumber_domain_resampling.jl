##### This is experiemental code for resampling in the frequency domain, not functional

## freq domain null resampling
function multitaper_estimate_resampled(
        data,
        region;
        nfreq,
        fmax,
        tapers,
        mean_method::MeanEstimationMethod = DefaultMean(),
        nresamples::Int
)
    data, dim = check_spatial_data(data)
    freq = make_freq(nfreq, fmax, dim)
    J_n = tapered_dft(data, tapers, nfreq, fmax, region, mean_method)
    power = dft2spectralmatrix(J_n)
    resampled = [SpectralEstimate(freq, dft2spectralmatrix(null_resample(J_n)))
                 for _ in 1:nresamples]
    observed = SpectralEstimate(freq, power)
    return (observed = observed, resampled = resampled)
end

"""
	null_resample(J_n::Array)

Resamples the DFTs across M independently in P, assuming that the DFTs are stored as one large array which is P x M x n_1 x ... x n_D.
"""
function null_resample(J_n::Array{T, N}) where {T, N}
    error("under development!!")
    M = size(J_n, 2)
    P = size(J_n, 1)
    permutedims(
        map(j -> selectdim(selectdim(J_n, 1, j), 2, rand(1:M, M)), 1:P),
        (N, 1:(N - 1)...)
    )
end

"""
	null_resample(J_n::NTuple{P, Array{T, N}}) where {P, T, N}

Resamples the DFTs across M independently in P, assuming that the DFTs are stored as a tuple of P arrays of size n_1 x ... x n_D x M.
"""
function null_resample(J_n::NTuple{P, Array{T, N}}) where {P, T, N}
    M = size(first(J_n))[end]
    # ntuple(j -> selectdim(J_n[j], N, rand(1:M, M)), Val{P}())
    new_J = ntuple(p -> Array{T, N}(undef, size(J_n[p])), Val{P}())
    for p in 1:P
        for i in CartesianIndices(J_n[p])
            new_J[p][i] = J_n[p][i.I[1:(end - 1)]..., rand(1:M)]
        end
    end
    return new_J
end
