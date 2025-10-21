struct PartialResampler{P, C, M, K, R}
    precompute::P
    cross_resampler::C
    marginal_resampler::M
    mt_kwargs::K
    region::R
end

"""
    PartialResampler(
        data::NTuple{P, PointSet},
        region;
        tapers,
        nk,
        kmax,
        mean_method = DefaultMean(),
        shift_method,
        nk_marginal_compute,
        kmax_marginal_compute
) where {P}

Creates a `PartialResampler` object that can be used to generate partial spectra.
"""
PartialResampler(data, region; kwargs...) = PartialResampler(
    spatial_data(data, region); kwargs...)
function PartialResampler(
        data::SpatialData;
        tapers,
        nk,
        kmax,
        mean_method = DefaultMean(),
        shift_method,
        nk_marginal_compute,
        kmax_marginal_compute
)
    precompute = create_resampler_precompute(
        data; nk = nk, kmax = kmax, tapers = tapers, mean_method = mean_method
    )
    cross_resampler = PartialCrossResampler(
        data, shift_method
    )
    marginal_resampler = PartialMarginalResampler(
        data; tapers = tapers, nk = nk_marginal_compute,
        kmax = kmax_marginal_compute
    )
    mt_kwargs = (tapers = tapers, nk = nk, kmax = kmax, mean_method = mean_method)
    return PartialResampler(
        precompute, cross_resampler, marginal_resampler, mt_kwargs, getregion(data))
end

function Base.rand(rng::AbstractRNG, resampler::PartialResampler)
    region = resampler.region
    tapers = resampler.mt_kwargs.tapers
    nk = resampler.mt_kwargs.nk
    kmax = resampler.mt_kwargs.kmax
    mean_method = resampler.mt_kwargs.mean_method

    d_cross = rand(rng, resampler.cross_resampler)
    d_marginal = rand(rng, resampler.marginal_resampler)

    J_cross = tapered_dft(d_cross, tapers, nk, kmax, mean_method[1:(end - 1)])
    mean_cross = mean_estimate(d_cross, mean_method[1:(end - 1)])
    J_marginal = tapered_dft(d_marginal, tapers, nk, kmax, mean_method)
    mean_marginal = mean_estimate(d_marginal, mean_method)
    atoms = covariance_zero_atom(d_marginal) # because cross doesn't have atoms

    return partial_from_resampled(
        J_cross, mean_cross,
        J_marginal, mean_marginal,
        atoms,
        resampler.precompute.wavenumber,
        resampler.precompute.J_original,
        resampler.precompute.mean_original,
        resampler.precompute.f_inv_cross,
        resampler.precompute.f_inv_marginal,
        resampler.precompute.ntapers
    )
end

function create_resampler_precompute(
        data::MultipleSpatialDataTuple; nk, kmax, tapers, mean_method = DefaultMean())
    _nk, _kmax = _validate_wavenumber_params(nk, kmax, nothing, data)
    wavenumber = _make_wavenumber_grid(_nk, _kmax)
    J_n = tapered_dft(data, tapers, _nk, _kmax, mean_method)
    power = _dft_to_spectral_matrix(data, J_n, _nk)
    f_inv_cross = create_f_inv_cross(power)
    f_inv_marginal = create_f_inv_marginal(power)
    mean_original = mean_estimate(data, mean_method)

    return (wavenumber = wavenumber, J_original = J_n, f_inv_cross = f_inv_cross,
        f_inv_marginal = f_inv_marginal,
        mean_original = mean_original, ntapers = length(tapers))
end

function create_f_inv_marginal(power::AbstractArray{
        SMatrix{P, P, T, S}, D}) where {P, T, S, D}
    return ntuple(
        p -> inv.(getindex.(
            power, Ref(static_not(Val{P}(), p)), Ref(static_not(Val{P}(), p)))),
        Val{P}()
    )
end

function create_f_inv_cross(power::AbstractArray{SMatrix{P, P, T, S}, D}) where {P, T, S, D}
    f_inv_cross = [zeros(SVector{P - 2, eltype(eltype(power))}, size(power))
                   for p in 1:(P - 1), q in 1:(P - 1)]

    for p in 1:(P - 1), q in (p + 1):P
        idx = static_not(Val{P}(), p, q)
        f_inv_cross[p, q - p] = getindex.(power, Ref(idx), Ref(idx)) .\
                                getindex.(power, Ref(idx), Ref(q)) # premultiplied form
    end
    return f_inv_cross
end

function create_f_inv_cross(power::AbstractArray{SMatrix{2, 2, T, S}, D}) where {T, S, D}
    return nothing
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
        wavenumber,
        J_original,
        mean_original,
        f_inv_cross,
        f_inv_marginal,
        ntapers::Int) where {Q}
        J_cross::NTuple{Q, AbstractArray{T, N}},
        J_marginal::NTuple{P, AbstractArray{T, N}},
        atoms,
        wavenumber,
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
        J_cross::AbstractArray{<:SVector{Q}},
        mean_cross,
        J_marginal,
        mean_marginal,
        atoms,
        wavenumber,
        J_original,
        mean_original,
        f_inv_cross,
        f_inv_marginal,
        ntapers::Int) where {Q}
    par_spec = partial_from_resampled(
        J_cross, mean_cross, J_marginal, mean_marginal, atoms, wavenumber, J_original,
        mean_original, f_inv_cross, f_inv_marginal, nothing
    )
    par = get_estimates(par_spec)
    for i in CartesianIndices(par)
        denom = ntapers .* ones(typeof(par[i])) .- Q .+ 2 - I
        par[i] = par[i] .* ntapers ./ denom
    end
    return Spectra{PartialTrait}(
        wavenumber, par, get_process_information(par_spec), EstimationInformation(ntapers))
end

function partial_from_resampled(
        J_cross::AbstractArray{<:SVector{Q}, N},
        mean_cross,
        J_marginal::AbstractArray{<:SVector{P}, N},
        mean_marginal,
        atoms,
        wavenumber,
        J_original::AbstractArray{<:SVector{P}, N},
        mean_original,
        f_inv_cross::Union{Matrix{<:AbstractArray{SVector{R, T}, D}}, Nothing},
        f_inv_marginal::NTuple{P, AbstractArray{SMatrix{Q, Q, T, L}, D}},
        ntapers::Nothing) where {P, Q, R, L, T, N, D}
    @assert Q == P - 1
    @assert (isnothing(f_inv_cross) && P == 2) || R == P - 2
    @assert isnothing(f_inv_cross) || size(f_inv_cross) == (P - 1, P - 1)

    output = zeros(SMatrix{P, P, T, P^2}, size(J_cross)[1:(end - 1)]) # last dimension is tapers
    for i in CartesianIndices(output) # todo: now very outdated
        for p in 1:P
            other_idx = static_not(Val{P}(), p)
            f_pp = mean(
                abs2, J_marginal[i, m][p] for m in axes(J_marginal, ndims(J_marginal)))
            f_px = mean(J_marginal[i, m][p] * J_original[i, m][other_idx]'
            for m in axes(J_marginal, N))
            x = f_pp - f_px * f_inv_marginal[p][i] * f_px'
            output[i] = setindex(output[i], x, p, p)
        end
        for p in 1:(P - 1), q in (p + 1):P
            x = partial_from_resampled_cross(J_cross, J_original, f_inv_cross, p, q, i)
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

    process_information = ProcessInformation{D, MultipleTupleTrait}(
        1:P, 1:P, mean_product, atoms) # should get the means for each part separately, means we should have means that are a matrix of products
    return Spectra{PartialTrait}(
        wavenumber, output, process_information, EstimationInformation(ntapers))
end

function partial_from_resampled_cross(
        J_cross, J_original::AbstractArray{<:SVector{P}}, f_inv_cross, p, q, i) where {P}
    other_idx = static_not(Val{P}(), p, q)
    f_pq = mean(J_cross[i, m][p] * conj(J_original[i, m][q])
    for m in axes(J_cross, ndims(J_cross)))
    f_px = mean(J_cross[i, m][p] * J_original[i, m][other_idx]'
    for m in axes(J_cross, ndims(J_cross)))
    return f_pq - f_px * f_inv_cross[p, q - p][i]
end

function partial_from_resampled_cross(
        J_cross, J_original::AbstractArray{<:SVector{2}}, f_inv_cross::Nothing, p, q, i)
    f_pq = mean(J_cross[i, m][p] * conj(J_original[i, m][q])
    for m in axes(J_cross, ndims(J_cross)))
    return f_pq
end

"""
    partial_shift_resample(rng::AbstractRNG, statistic, resampler; kwargs...)

Assumes that `statistic` is of the form
    statistic(partialspectra::PartialSpectra; kwargs...)

The `resampler` should be able to generate resamples of the partial spectra.
"""
function partial_shift_resample(rng::AbstractRNG, statistic, resampler; kwargs...)
    partial_spec = rand(rng, resampler)
    return statistic(partial_spec; kwargs...)
end
function partial_shift_resample(statistic, resampler; kwargs...)
    partial_shift_resample(Random.default_rng(), statistic, resampler; kwargs...)
end
# TODO: currently we rotationally average K function in usual but not in the resampling, should think about the best way to handle this
