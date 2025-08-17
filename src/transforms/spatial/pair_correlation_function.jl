struct PairCorrelationFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    paircorrelation_function::T
    function PairCorrelationFunction(radii::R, pcf::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, pcf)
        new{R,T,D,P}(radii, pcf)
    end
end

getargument(f::PairCorrelationFunction) = f.radii
getestimate(f::PairCorrelationFunction) = f.paircorrelation_function
getextrafields(::PairCorrelationFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

abstract type PCFMethod end
struct PCFMethodA <: PCFMethod end
struct PCFMethodB <: PCFMethod end
struct PCFMethodC <: PCFMethod end
struct PCFMethodD <: PCFMethod end

function K2paircorrelation(radii, k, ::Val{D}, penalty, ::PCFMethodA) where {D}
    A = 2 * (π)^(D / 2) / gamma(D/2)
    kinterp = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k, penalty)
    ∂k = BSplineKit.Derivative(1) * kinterp
    return ∂k.(radii) ./ (A.*radii.^(D-1))
end

function K2paircorrelation(radii, k, ::Val{2}, penalty, ::PCFMethodA)
    kinterp = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k, penalty)
    ∂k = BSplineKit.Derivative(1) * kinterp
    return ∂k.(radii) ./ (2pi.*radii)
end

function K2paircorrelation(radii, k, ::Val{2}, penalty, ::PCFMethodB)
    y = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k ./ (2π.*radii), penalty)
    ∂y = BSplineKit.Derivative(1) * y
    return ∂y.(radii) + k ./ (2π.*radii.^2)
end

function K2paircorrelation(radii, k, ::Val{2}, penalty, ::PCFMethodC)
    z = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, k ./ (2π.*radii.^2), penalty)
    ∂z = BSplineKit.Derivative(1) * z
    return ∂z.(radii) .* radii + k ./ (π.*radii.^2)
end

function K2paircorrelation(radii, k, ::Val{2}, penalty, ::PCFMethodD)
    v = BSplineKit.fit(BSplineKit.BSplineOrder(4), radii, sqrt.(k), penalty)
    ∂v = BSplineKit.Derivative(1) * v
    return ∂v.(radii) .* sqrt.(k) ./ (π.*radii)
end

function paircorrelation_function(k::KFunction{R,T,D}; penalty = 0.0, method = PCFMethodC()) where {R,T,D}
    pcf = Dict(index => K2paircorrelation(k.radii, val, Val{D}(), penalty, method) for (index, val) in k.K_function)
    return PairCorrelationFunction(k.radii, pcf, Val{D}())
end

function paircorrelation_function(
    data,
    region,
    radii,
    indices = default_indices(data);
    penalty = 0.0,
    method::PCFMethod = PCFMethodC(),
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    k = K_function(
        data,
        region,
        radii,
        indices;
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
    )
    return paircorrelation_function(k; penalty = penalty, method = method)
end

## direct method

function paircorrelation_function_direct(
    f::SpectralEstimate{D,F,P,N},
    λ,
    radii,
    indices = default_indices(f),
) where {D,F,P,N}
    check_spectral_indices(indices, f)
    pcf = Dict(
        index => sdf2pcf(f[index...], λ[index[1]], λ[index[2]], radii) for index in indices
    )
    return PairCorrelationFunction(radii, pcf, Val{D}())
end

function paircorrelation_function_direct(
    data,
    region,
    radii,
    indices = default_indices(data);
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    check_spectral_indices(indices, data)
    f = multitaper_estimate(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        mean_method = mean_method,
    )
    λ = mean_estimate(data, region, mean_method)
    return paircorrelation_function_direct(f, λ, radii, indices)
end


function sdf2pcf(f, λ1, λ2, radii::AbstractVector{<:Number})
    [sdf2pcf(f, λ1, λ2, radius) for radius in radii]
end

"""
    sdf2pcf(f, λ1, λ2, radius::Number)

Takes some form of spectra and returns the C function for the `radius`.
"""
function sdf2pcf(
    f::Union{SpectralEstimate{D,F,P,N},PartialSpectra{D,F,P,N}},
    λ1,
    λ2,
    radius::Number,
) where {D,F,P,N}
    freq = getargument(f)
    spectra = getestimate(f)
    pcf_unweighted =
        prod(step, freq) * real(
            sum(
                f * pcf_weight(radius, k, Val{D}()) for
                (f, k) in zip(spectra, Iterators.product(freq...))
            ),
        )
    return pcf_unweighted ./ (λ1 * λ2) .+ 1
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