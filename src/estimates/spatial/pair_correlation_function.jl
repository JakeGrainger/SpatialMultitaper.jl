struct PairCorrelationFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function PairCorrelationFunction{E}(
            radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        P, Q = checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        new{E, D, P, Q, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end
getargument(f::PairCorrelationFunction) = f.radii
getestimate(f::PairCorrelationFunction) = f.value

function pair_correlation_function(data, region; kwargs...)
    pair_correlation_function(spatial_data(data, region); kwargs...)
end
function pair_correlation_function(data::SpatialData; pcf_method = PCFMethodC(), kwargs...)
    pair_correlation_function(k_function(data; kwargs...); pcf_method = pcf_method)
end
function pair_correlation_function(est::KFunction{E}; pcf_method = PCFMethodC()) where {E}
    radii = getargument(est)
    value = k2paircorrelation(est, pcf_method)
    processinfo = getprocessinformation(est)
    estimationinfo = getestimationinformation(est)
    return PairCorrelationFunction{E}(radii, value, processinfo, estimationinfo)
end
function pair_correlation_function(spectrum::Spectra; kwargs...)
    pair_correlation_function(k_function(spectrum); kwargs...)
end

function partial_pair_correlation_function(data, region; kwargs...)
    partial_pair_correlation_function(spatial_data(data, region); kwargs...)
end
function partial_pair_correlation_function(
        data::SpatialData; pcf_method = PCFMethodC(), kwargs...)
    pair_correlation_function(
        partial_k_function(data; kwargs...); pcf_method = pcf_method)
end
function partial_pair_correlation_function(spectrum::Spectra{MarginalTrait}; kwargs...)
    pair_correlation_function(partial_spectra(spectrum); kwargs...)
end
function partial_pair_correlation_function(spectrum::Spectra{PartialTrait}; kwargs...)
    pair_correlation_function(spectrum; kwargs...)
end
function partial_pair_correlation_function(est::KFunction{PartialTrait}; kwargs...)
    pair_correlation_function(est; kwargs...)
end
function partial_pair_correlation_function(est::CFunction{PartialTrait}; kwargs...)
    throw(partial_from_marginal_error(PairCorrelationFunction, typeof(est)))
end

## direct method
function pair_correlation_function_direct(data, region; kwargs...)
    pair_correlation_function_direct(spatial_data(data, region); kwargs...)
end
function pair_correlation_function_direct(data::SpatialData; radii, spectra_kwargs...)
    spectrum = spectra(data; spectra_kwargs...)
    return pair_correlation_function_direct(spectrum, radii = radii)
end

function pair_correlation_function_direct(f::Spectra{E}; radii) where {E}
    value = sdf2pcf(f, radii)
    processinfo = getprocessinformation(f)
    estimationinfo = getestimationinformation(f)
    return PairCorrelationFunction{E}(radii, value, processinfo, estimationinfo)
end

function partial_pair_correlation_function_direct(data, region; kwargs...)
    partial_pair_correlation_function_direct(spatial_data(data, region); kwargs...)
end
function partial_pair_correlation_function_direct(
        data::SpatialData; radii, spectra_kwargs...)
    spectrum = partial_spectra(data; spectra_kwargs...)
    return pair_correlation_function_direct(spectrum, radii = radii)
end
function partial_pair_correlation_function_direct(spectrum::Spectra{PartialTrait}; radii)
    return pair_correlation_function_direct(spectrum, radii = radii)
end

function partial_pair_correlation_function_direct(spectrum::Spectra{MarginalTrait}; radii)
    return pair_correlation_function_direct(partial_spectra(spectrum), radii = radii)
end

## internals

abstract type PCFMethod end
@kwdef struct PCFMethodA{T} <: PCFMethod
    penalty::T = 0.0
end
@kwdef struct PCFMethodB{T} <: PCFMethod
    penalty::T = 0.0
end
@kwdef struct PCFMethodC{T} <: PCFMethod
    penalty::T = 0.0
end
@kwdef struct PCFMethodD{T} <: PCFMethod
    penalty::T = 0.0
end

function k2paircorrelation(est::KFunction{E, D, P, Q}, method) where {E, D, P, Q}
    radii = getargument(est)
    stacked_pcf = stack(_k2paircorrelation(radii, getestimate(est[i, j]), Val{D}(), method)
    for i in 1:P, j in 1:Q)
    return _post_process_pcf(stacked_pcf, est)
end

function _post_process_pcf(stacked_pcf, template_est::KFunction)
    return _reshape_pcf_output(stacked_pcf, getestimate(template_est), P, Q)
end

function _reshape_pcf_output(data, ::AbstractVector{<:SMatrix{P, Q}}) where {P, Q}
    [SMatrix{P, Q, T, L}(data[:, :, k]) for k in axes(data, 3)]
end

_reshape_pcf_output(data, ::AbstractArray{<:Number, 3}, P, Q) = permutedims(data, (2, 3, 1))

function _k2paircorrelation(radii, k, ::Val{D}, method::PCFMethodA) where {D}
    penalty = method.penalty
    A = 2 * (π)^(D / 2) / gamma(D / 2)
    kinterp = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k, penalty)
    ∂k = BSplineKit.Derivative(1) * kinterp
    return ∂k.(radii) ./ (A .* radii .^ (D - 1))
end

function _k2paircorrelation(radii, k, ::Val{2}, method::PCFMethodA)
    penalty = method.penalty
    kinterp = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k, penalty)
    ∂k = BSplineKit.Derivative(1) * kinterp
    return ∂k.(radii) ./ (2pi .* radii)
end

function _k2paircorrelation(radii, k, ::Val{2}, method::PCFMethodB)
    penalty = method.penalty
    y = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k ./ (2π .* radii), penalty)
    ∂y = BSplineKit.Derivative(1) * y
    return ∂y.(radii) + k ./ (2π .* radii .^ 2)
end

function _k2paircorrelation(radii, k, ::Val{2}, method::PCFMethodC)
    penalty = method.penalty
    z = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k ./ (2π .* radii .^ 2), penalty)
    ∂z = BSplineKit.Derivative(1) * z
    return ∂z.(radii) .* radii + k ./ (π .* radii .^ 2)
end

function _k2paircorrelation(radii, k, ::Val{2}, method::PCFMethodD)
    penalty = method.penalty
    v = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, sqrt.(k), penalty)
    ∂v = BSplineKit.Derivative(1) * v
    return ∂v.(radii) .* sqrt.(k) ./ (π .* radii)
end

## dirtect method internals

function sdf2pcf(f, radii::AbstractVector{<:Number})
    [_sdf2pcf(f, radius) for radius in radii]
end

function _sdf2pcf(
        f::Spectra{E, D},
        radius::Number
) where {E, D}
    freq = getargument(f)
    spectra = getestimate(f)
    zeroatom = getprocessinformation(f).atoms
    mean_prod = getprocessinformation(f).mean_product
    pcf_unweighted = prod(step, freq) * real(
        sum((f - zeroatom) * pcf_weight(radius, k, Val{D}())
    for
    (f, k) in zip(spectra, Iterators.product(freq...))
    ))
    return pcf_unweighted ./ (mean_prod) .+ 1
end

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
