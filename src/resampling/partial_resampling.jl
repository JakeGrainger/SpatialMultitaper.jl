struct PartialResampler{P, C, M}
    precompute::P
    cross_resampler::C
    marginal_resampler::M
end

function Base.rand(rng::AbstractRNG, resampler::PartialResampler)
    partial_from_resampled(
        rand(rng, resampler.cross_resampler),
        rand(rng, resampler.marginal_resampler),
        resampler.precompute.freq,
        resampler.precompute.J_original,
        resampler.precompute.f_inv_cross,
        resampler.precompute.f_inv_marginal,
        resampler.precompute.ntapers
    )
end

function create_resampler_precompute(
        data::NTuple{P}, region; nfreq, fmax, tapers, mean_method = DefaultMean()) where {P}
    mask.(data, Ref(region))
    data, dim = check_spatial_data(data)
    mean_method = check_mean_method(mean_method, data)
    J_n = tapered_dft(data, tapers, nfreq, fmax, region, mean_method)
    freq = make_freq(nfreq, fmax, dim)
    power = dft2spectralmatrix(J_n)

    f_inv_cross = ones(
        Array{SVector{P - 2, eltype(power)}, ndims(power)}, P - 1, P - 1)
    for p in 1:(P - 1), q in (p + 1):P
        idx = static_not(Val{P}(), p, q)
        f_inv_cross[p, q - p] = getindex.(power, Ref(idx), Ref(idx)) .\
                                getindex.(power, Ref(idx), Ref(q)) # premultiplied form
    end

    f_inv_marginal = ntuple(
        p -> inv.(getindex.(
            power, Ref(static_not(Val{P}(), p)), Ref(static_not(Val{P}(), p)))),
        Val{P - 1}()
    )

    return (freq = freq, J_original = J_n, f_inv_cross = f_inv_cross,
        f_inv_marginal = f_inv_marginal, ntapers = length(tapers))
end

function static_not(::Val{P}, p) where {P}
    StaticArrays.sacollect(
        SVector{P - 1, Int}, ApplyArray(vcat, 1:(p - 1), (p + 1):P))
end
function static_not(::Val{P}, p, q) where {P}
    @assert p < q
    StaticArrays.sacollect(
        SVector{P - 2, Int}, ApplyArray(vcat, 1:(p - 1), (p + 1):(q - 1), (q + 1):P))
end

"""
    partial_from_resampled(
        J_cross::NTuple{Q, AbstractArray{T, N}},
        J_marginal::NTuple{P, AbstractArray{T, N}},
        freq,
        J_original::NTuple{P, AbstractArray{T, N}},
        f_inv_cross::Matrix{AbstractArray{SMatrix{R, Q, T, S}, D}}},
        f_inv_marginal::NTuple{P, AbstractArray{SMatrix{Q, Q, T, L}, D}},
        ntapers::Int) where {P, Q, R, L, S, T, N, D}

Takes the original DFTs, the inverses of the cross and marginal spectral estimates from the
original data, and the DFTs from the resampled data, and returns the partial spectra.

The cross inverses will have been premultiplied so that they take the form
f[idx, idx]^{-1}f[idx, jdx]
"""
function partial_from_resampled(
        J_cross::NTuple{Q},
        J_marginal,
        freq,
        J_original,
        f_inv_cross,
        f_inv_marginal,
        ntapers::Int) where {Q}
    par_spec = partial_from_resampled(
        J_cross, J_marginal, freq, J_original, f_inv_cross, f_inv_marginal, nothing
    )
    par = par_spec.partial_spectra
    for i in CartesianIndices(par)
        denom = ntapers .* ones(typeof(par[i])) .- Q .+ 2 - I
        par[i] = par[i] .* ntapers ./ denom
    end
    return PartialSpectra(freq, par)
end

function partial_from_resampled(
        J_cross::NTuple{Q, AbstractArray{T, N}},
        J_marginal::NTuple{P, AbstractArray{T, N}},
        freq,
        J_original::NTuple{P, AbstractArray{T, N}},
        f_inv_cross::NTuple{P, AbstractArray{SMatrix{R, Q, T, S}, D}},
        f_inv_marginal::NTuple{P, AbstractArray{SMatrix{Q, Q, T, L}, D}},
        ntapers::Nothing) where {P, Q, R, L, S, T, N, D}
    @assert Q == P - 1
    @assert R == P - 2

    output = zeros(SMatrix{P, P, T, P^2}, size(J_cross[1]))
    for i in CartesianIndices(S_mat)
        for p in 1:P
            other_idx = ApplyArray(vcat, 1:(p - 1), (p + 1):P)
            J_p = J_marginal[p]
            f_pp = mean(abs2, view(J_p, i, :))
            f_px = mean(J_p[i, m] *
                        SVector(ntuple(j -> J_original[other_idx[j]][i, m], Val{P - 1}()))'
            for m in axes(J_p, N))
            x = f_pp - f_px * f_inv_marginal[p][i] * f_px'
            setindex(S_mat[i], x, p, p)
        end
        for p in 1:(P - 1), q in (p + 1):P
            other_idx = ApplyArray(vcat, 1:(p - 1), (p + 1):(q - 1), (q + 1):P)
            J_p = J_marginal[p]
            J_q = J_marginal[q]
            f_pq = mean(J_p[i, m] * conj(J_q[i, m]) for m in axes(J_p, N))
            f_px = mean(J_p[i, m] *
                        SVector(ntuple(j -> J_original[other_idx[j]][i, m], Val{P - 2}()))'
            for m in axes(J_p, N))
            x = f_pq - f_px * f_inv_cross[p, q - p][i]
            setindex(S_mat[i], x, p, q)
            setindex(S_mat[i], conj(x), q, p)
        end
    end
    return PartialSpectra(freq, output)
end
