struct PairCorrelationFunction{R,P}
    radii::R
    pcf::P
end

function paircorrelation_function(
    f::SpectralEstimate,
    λ,
    radii,
    indices = default_indices(f),
)
    check_spectral_indices(indices, f)
    pcf = Dict(
        index => sdf2pcf(f[index...], λ[index[1]], λ[index[2]], radii) for index in indices
    )
    return PairCorrelationFunction(radii, pcf)
end

function paircorrelation_function(
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
    λ = mean_estimate.(data, Ref(region), Ref(mean_method)) # TODO: technically this isn't compatible with mean methods as they should be passed as Tuples.
    return paircorrelation_function(f, λ, radii, indices)
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
    freq = getfreq(f)
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