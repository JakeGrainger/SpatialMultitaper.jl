struct PartialResampler{P, C, M, K, R}
    precompute::P
    cross_resampler::C
    marginal_resampler::M
    mt_kwargs::K
    region::R
end

function PartialResampler(
        data::NTuple{P, PointSet},
        region;
        tapers,
        nfreq,
        fmax,
        mean_method = DefaultMean(),
        shift_method,
        nfreq_marginal_compute,
        fmax_marginal_compute
) where {P}
    precompute = create_resampler_precompute(
        data, region; nfreq = nfreq, fmax = fmax, tapers = tapers, mean_method = mean_method
    )
    cross_resampler = PartialCrossResampler(
        data, shift_method
    )
    marginal_resampler = PartialMarginalResampler(
        data, region; tapers = tapers, nfreq = nfreq_marginal_compute,
        fmax = fmax_marginal_compute
    )
    mt_kwargs = (tapers = tapers, nfreq = nfreq, fmax = fmax, mean_method = mean_method)
    return PartialResampler(
        precompute, cross_resampler, marginal_resampler, mt_kwargs, region)
end

function Base.rand(rng::AbstractRNG, resampler::PartialResampler)
    region = resampler.region
    tapers = resampler.mt_kwargs.tapers
    nfreq = resampler.mt_kwargs.nfreq
    fmax = resampler.mt_kwargs.fmax
    mean_method = resampler.mt_kwargs.mean_method

    d_cross = rand(rng, resampler.cross_resampler)
    d_marginal = rand(rng, resampler.marginal_resampler)

    J_cross = tapered_dft(d_cross, tapers, nfreq, fmax, region, mean_method[1:(end - 1)])
    mean_cross = mean_estimate(d_cross, region, mean_method[1:(end - 1)])
    J_marginal = tapered_dft(d_marginal, tapers, nfreq, fmax, region, mean_method)
    mean_marginal = mean_estimate(d_marginal, region, mean_method)
    atoms = atom_estimate(d_marginal, region) # because cross doesn't have atoms

    return partial_from_resampled(
        J_cross, mean_cross,
        J_marginal, mean_marginal,
        atoms,
        resampler.precompute.freq,
        resampler.precompute.J_original,
        resampler.precompute.mean_original,
        resampler.precompute.f_inv_cross,
        resampler.precompute.f_inv_marginal,
        resampler.precompute.ntapers
    )
end

function create_resampler_precompute(
        data::NTuple{P}, region; nfreq, fmax, tapers, mean_method = DefaultMean()) where {P}
    mask.(data, Ref(region))
    data, dim = check_spatial_data(data)
    J_n = tapered_dft(data, tapers, nfreq, fmax, region, mean_method)
    freq = make_freq(nfreq, fmax, dim)
    power = dft2spectralmatrix(J_n)

    f_inv_cross = [zeros(SVector{P - 2, eltype(eltype(power))}, size(power))
                   for p in 1:(P - 1), q in 1:(P - 1)]
    for p in 1:(P - 1), q in (p + 1):P
        idx = static_not(Val{P}(), p, q)
        f_inv_cross[p, q - p] = getindex.(power, Ref(idx), Ref(idx)) .\
                                getindex.(power, Ref(idx), Ref(q)) # premultiplied form
    end

    f_inv_marginal = ntuple(
        p -> inv.(getindex.(
            power, Ref(static_not(Val{P}(), p)), Ref(static_not(Val{P}(), p)))),
        Val{P}()
    )

    mean_original = mean_estimate(data, region, mean_method)

    return (freq = freq, J_original = J_n, f_inv_cross = f_inv_cross,
        f_inv_marginal = f_inv_marginal, mean_original = mean_original, ntapers = length(tapers))
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
        J_cross::NTuple{Q},
        mean_cross,
        J_marginal,
        mean_marginal,
        atoms,
        freq,
        J_original,
        mean_original,
        f_inv_cross,
        f_inv_marginal,
        ntapers::Int) where {Q}
        J_cross::NTuple{Q, AbstractArray{T, N}},
        J_marginal::NTuple{P, AbstractArray{T, N}},
        atoms,
        freq,
        J_original::NTuple{P, AbstractArray{T, N}},
        f_inv_cross::Matrix{AbstractArray{SMatrix{R, Q, T, S}, D}}},
        f_inv_marginal::NTuple{P, AbstractArray{SMatrix{Q, Q, T, L}, D}},
        λ,
        ntapers::Int) where {P, Q, R, L, S, T, N, D}

Takes the original DFTs, the inverses of the cross and marginal spectral estimates from the
original data, and the DFTs from the resampled data, and returns the partial spectra.

The cross inverses will have been premultiplied so that they take the form
f[idx, idx]^{-1}f[idx, jdx]
"""
function partial_from_resampled(
        J_cross::NTuple{Q},
        mean_cross,
        J_marginal,
        mean_marginal,
        atoms,
        freq,
        J_original,
        mean_original,
        f_inv_cross,
        f_inv_marginal,
        ntapers::Int) where {Q}
    par_spec = partial_from_resampled(
        J_cross, mean_cross, J_marginal, mean_marginal, atoms, freq, J_original,
        mean_original, f_inv_cross, f_inv_marginal, nothing
    )
    par = par_spec.partial_spectra
    for i in CartesianIndices(par)
        denom = ntapers .* ones(typeof(par[i])) .- Q .+ 2 - I
        par[i] = par[i] .* ntapers ./ denom
    end
    return PartialSpectra(
        freq, par, getprocessinformation(par_spec), EstimationInformation(ntapers))
end

function partial_from_resampled(
        J_cross::NTuple{Q, AbstractArray{T, N}},
        mean_cross,
        J_marginal::NTuple{P, AbstractArray{T, N}},
        mean_marginal,
        atoms,
        freq,
        J_original::NTuple{P, AbstractArray{T, N}},
        mean_original,
        f_inv_cross::Matrix{<:AbstractArray{SVector{R, T}, D}},
        f_inv_marginal::NTuple{P, AbstractArray{SMatrix{Q, Q, T, L}, D}},
        ntapers::Nothing) where {P, Q, R, L, T, N, D}
    @assert Q == P - 1
    @assert R == P - 2
    @assert size(f_inv_cross) == (P - 1, P - 1)

    output = zeros(SMatrix{P, P, T, P^2}, size(J_cross[1])[1:(end - 1)]) # last dimension is tapers
    for i in CartesianIndices(output)
        for p in 1:P
            other_idx = ApplyArray(vcat, 1:(p - 1), (p + 1):P)
            J_p = J_marginal[p]
            f_pp = mean(abs2, view(J_p, i, :))
            f_px = mean(J_p[i, m] *
                        SVector(ntuple(j -> J_original[other_idx[j]][i, m], Val{P - 1}()))'
            for m in axes(J_p, N))
            x = f_pp - f_px * f_inv_marginal[p][i] * f_px'
            output[i] = setindex(output[i], x, p, p)
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
            output[i] = setindex(output[i], x, p, q)
            output[i] = setindex(output[i], conj(x), q, p)
        end
    end
    mean_product = zero(SMatrix{P, P, eltype(mean_original), P^2})
    for p in 1:P
        mean_product = setindex(mean_product, mean_marginal[p]^2, p, p)
    end
    for p in 1:(P - 1), q in (p + 1):P
        x = mean_cross[p] * mean_original[q]
        mean_product = setindex(mean_product, x, p, q)
        mean_product = setindex(mean_product, x, q, p)
    end

    process_information = ProcessInformation(1:P, 1:P, mean_product, atoms, Val{D}()) # should get the means for each part separately, means we should have means that are a matrix of products
    return PartialSpectra(
        freq, output, process_information, EstimationInformation(ntapers))
end
