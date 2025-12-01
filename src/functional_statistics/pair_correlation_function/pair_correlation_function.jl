struct PairCorrelationFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function PairCorrelationFunction{E}(
            radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        new{E, D, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end

## required interface
function computed_from(::Type{<:PairCorrelationFunction{E, D}}) where {E, D}
    (KFunction{E, D}, Spectra{E, D})
end

function allocate_estimate_memory(::Type{<:PairCorrelationFunction}, ::Type{S},
        relevant_memory; kwargs...) where {S <: Spectra}
    mem = relevant_memory[1:2]
    wavenumber = _extract_wavenumber_from_c_mem(S, relevant_memory[3]; kwargs...)
    spatial_output = preallocate_radial_output(mem...; kwargs...)
    weights = precompute_pcf_weights(S, mem[1], wavenumber; kwargs...)
    return spatial_output, weights
end
function allocate_estimate_memory(
        ::Type{<:PairCorrelationFunction}, ::Type{S},
        relevant_memory; kwargs...) where {S <: KFunction}
    return relevant_memory[1], nothing
end

function extract_relevant_memory(
        ::Type{<:PairCorrelationFunction}, est::AbstractEstimate)
    return deepcopy(get_estimates(est)), process_trait(est), get_evaluation_points(est)
end
function extract_relevant_memory(
        ::Type{<:PairCorrelationFunction}, mem::EstimateMemory)
    return mem.output_memory, process_trait(mem), nothing
end

function validate_core_parameters(::Type{<:PairCorrelationFunction}; kwargs...)
    if :radii in keys(kwargs)
        radii = kwargs[:radii]
        validate_radii(radii)
    end
    return nothing
end

function resolve_missing_parameters(
        ::Type{<:PairCorrelationFunction}, ::Type{<:Spectra}, arg; kwargs...)
    radii = get(kwargs, :radii, nothing)
    other_kwargs = filter(kv -> kv[1] != :radii, kwargs)
    return (radii = process_pcf_radii(radii, arg), other_kwargs...)
end

function resolve_missing_parameters(
        ::Type{<:PairCorrelationFunction}, ::Type{<:KFunction}, arg; kwargs...)
    pcf_method = get(kwargs, :pcf_method, nothing)
    radii = get(kwargs, :radii, nothing)
    other_kwargs = filter(kv -> kv[1] ∉ (:pcf_method, :radii), kwargs)
    return (pcf_method = process_pcf_method(pcf_method),
        radii = process_pcf_radii(radii, arg), other_kwargs...)
end

function validate_memory_compatibility(
        ::Type{<:PairCorrelationFunction}, mem, arg::Spectra; radii, kwargs...)
    validate_pcf_internal(mem.internal_memory, get_estimates(arg), process_trait(arg))
    validate_radial_memory(mem.output_memory, process_trait(arg), radii)
    return nothing
end

function validate_memory_compatibility(
        ::Type{<:PairCorrelationFunction}, mem, arg::KFunction; kwargs...)
    @argcheck size(mem.output_memory) == size(get_estimates(arg))
    @argcheck eltype(mem.output_memory) == eltype(get_estimates(arg))
    return nothing
end

function compute_estimate!(
        ::Type{<:PairCorrelationFunction{E}}, mem, arg::Spectra; radii, kwargs...) where {E}
    estimate = _sdf2pcf!(mem.output_memory, mem.internal_memory, arg, radii)

    processinfo = get_process_information(arg)
    estimationinfo = get_estimation_information(arg)
    return PairCorrelationFunction{E}(radii, estimate, processinfo, estimationinfo)
end

function compute_estimate!(::Type{<:PairCorrelationFunction{E}}, mem,
        arg::KFunction; pcf_method, kwargs...) where {E}
    radii = get_evaluation_points(arg)
    estimate = k2paircorrelation(arg, process_trait(arg), pcf_method) # TODO: modify to use memory

    processinfo = get_process_information(arg)
    estimationinfo = get_estimation_information(arg)
    return PairCorrelationFunction{E}(radii, estimate, processinfo, estimationinfo)
end

get_evaluation_points(f::PairCorrelationFunction) = f.radii

get_estimates(f::PairCorrelationFunction) = f.value

## additional interface

get_short_base_estimate_name(::Type{<:PairCorrelationFunction}) = "pcf"

## internals

### Validation

function validate_pcf_internal(weights, power, ::SingleProcessTrait)
    @argcheck size(weights)[1:(end - 1)] == size(power)
    return nothing
end

function validate_pcf_internal(weights, power, ::MultipleTupleTrait)
    @argcheck size(weights)[1:(end - 1)] == size(power)
    return nothing
end

function validate_pcf_internal(weights, power, ::MultipleVectorTrait)
    @argcheck size(weights)[1:(end - 1)] == size(power)[3:end]
    return nothing
end

### from spectra

function _sdf2pcf!(out, store, spectrum::Spectra, radii)
    power = get_estimates(spectrum)
    zero_atom = get_process_information(spectrum).atoms
    mean_prod = get_process_information(spectrum).mean_product
    trait = process_trait(spectrum)
    return _sdf2pcf_internal!(out, store, power, trait, zero_atom, mean_prod, radii)
end

function _sdf2pcf_internal!(
        out, store, power, ::Union{SingleProcessTrait, MultipleTupleTrait},
        zero_atom, mean_prod, radii::AbstractVector)
    for i in eachindex(out)
        out[i] = _sdf2C_sum(selectdim(store, ndims(store), i), power, zero_atom) ./
                 (mean_prod) .+ 1
    end
    return out
end

function _sdf2pcf_internal!(
        out, store, power::AbstractArray{<:Number, N}, ::MultipleVectorTrait,
        zero_atom, mean_prod, radii::AbstractVector) where {N}
    for i in axes(out, ndims(out))
        store_slice = selectdim(store, ndims(store), i)
        for idx in CartesianIndices(size(out)[1:2])
            power_slice = view(power, idx, ntuple(Returns(:), N - 2)...)
            zero_atom_slice = _slice_zero_atom(zero_atom, idx)
            out[idx, i] = _sdf2C_sum(store_slice, power_slice, zero_atom_slice) /
                          mean_prod[idx] + 1
        end
    end
    return out
end

function precompute_pcf_weights(::Type{<:Spectra{E, D}}, power, wavenumber;
        radii::AbstractVector, kwargs...) where {E, D}
    wavenumber_size = length.(wavenumber)
    wavenumber_spacing = prod(step, wavenumber)
    T = promote_type(float(eltype(radii)), float(typeof(wavenumber_spacing)))
    store = zeros(T, wavenumber_size..., length(radii))
    for (i, radius) in enumerate(radii)
        for (idx, k) in zip(
            CartesianIndices(wavenumber_size), Iterators.product(wavenumber...))
            store[idx, i] = pcf_weight(radius, k, Val{D}()) * wavenumber_spacing
        end
    end
    return store
end

precompute_pcf_weights(::Type{<:KFunction}, power, wavenumber; kwargs...) = nothing # KFunction-based PCF does not need weights

function pcf_weight(r, u, ::Val{1})
    error("Pair correlation function is not implemented 1D yet.")
end

function pcf_weight(r, u, ::Val{2})
    rx = r * norm(u)
    return (rx < 1e-10) ? 1.0 : besselj0(2π * rx)
end

function pcf_weight(r, u, ::Val{D}) where {D}
    rx = r * norm(u)
    return (rx < 1e-10) ? 1.0 :
           (gamma(D / 2) / (2π)^(D / 2)) * (1 / rx)^(D / 2) * besselj(D / 2 - 1, 2π * rx)
end

### from k function

abstract type PCFMethod end
@kwdef struct PCFMethodA{T} <: PCFMethod
    penalty::T = 0.5
end
@kwdef struct PCFMethodB{T} <: PCFMethod
    penalty::T = 0.5
end
@kwdef struct PCFMethodC{T} <: PCFMethod
    penalty::T = 0.5
end
@kwdef struct PCFMethodD{T} <: PCFMethod
    penalty::T = 0.5
end

process_pcf_method(::Nothing) = PCFMethodC()
process_pcf_method(method::PCFMethod) = method
function process_pcf_method(pcf_method)
    throw(ArgumentError("Invalid pcf_method provided, should be a `PCFMethod`, but $pcf_method provided."))
end

function k2paircorrelation(est::KFunction{E, D}, ::SingleProcessTrait, method) where {E, D}
    radii = get_evaluation_points(est)
    return _k2paircorrelation(radii, get_estimates(est), Val{D}(), method)
end

function k2paircorrelation(est::KFunction{E, D}, trait, method) where {E, D}
    radii = get_evaluation_points(est)
    stacked_pcf = stack(_k2paircorrelation(
                            radii, get_estimates(est[i, j]), Val{D}(), method)
    for i in 1:size(est, 1), j in 1:size(est, 2))
    return _post_process_pcf(stacked_pcf, est)
end

function _post_process_pcf(stacked_pcf, template_est::KFunction)
    return _reshape_pcf_output(stacked_pcf, get_estimates(template_est))
end

function _reshape_pcf_output(data, ::AbstractVector{<:SMatrix{P, Q}}) where {P, Q}
    [SMatrix{P, Q, eltype(data), P * Q}(data[k, :, :]) for k in axes(data, 1)]
end

_reshape_pcf_output(data, ::AbstractArray{<:Number, 3}) = permutedims(data, (2, 3, 1))

function _k2paircorrelation(radii, k, ::Val{D}, method::PCFMethodA) where {D}
    penalty = spar_to_lambda(BSplineKit.BSplineOrder(4), radii, k, method.penalty)
    kinterp = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k, penalty)
    A = 2 * (π)^(D / 2) / gamma(D / 2)
    ∂k = BSplineKit.Derivative(1) * kinterp
    return ∂k.(radii) ./ (A .* radii .^ (D - 1))
end

function _k2paircorrelation(radii, k, ::Val{2}, method::PCFMethodA)
    penalty = spar_to_lambda(BSplineKit.BSplineOrder(4), radii, k, method.penalty)
    kinterp = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k, penalty)
    ∂k = BSplineKit.Derivative(1) * kinterp
    return ∂k.(radii) ./ (2pi .* radii)
end

function _k2paircorrelation(radii, k, ::Val{2}, method::PCFMethodB)
    penalty = spar_to_lambda(
        BSplineKit.BSplineOrder(4), radii, k ./ (2π .* radii), method.penalty)
    y = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k ./ (2π .* radii), penalty)
    ∂y = BSplineKit.Derivative(1) * y
    return ∂y.(radii) + k ./ (2π .* radii .^ 2)
end

function _k2paircorrelation(radii, k, ::Val{2}, method::PCFMethodC)
    penalty = spar_to_lambda(
        BSplineKit.BSplineOrder(4), radii, k ./ (2π .* radii .^ 2), method.penalty)
    z = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k ./ (2π .* radii .^ 2), penalty)
    ∂z = BSplineKit.Derivative(1) * z
    return ∂z.(radii) .* radii + k ./ (π .* radii .^ 2)
end

function _k2paircorrelation(radii, k, ::Val{2}, method::PCFMethodD)
    penalty = spar_to_lambda(BSplineKit.BSplineOrder(4), radii, sqrt.(k), method.penalty)
    v = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, sqrt.(k), penalty)
    ∂v = BSplineKit.Derivative(1) * v
    return ∂v.(radii) .* sqrt.(k) ./ (π .* radii)
end
