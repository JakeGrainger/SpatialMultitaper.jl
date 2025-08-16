struct PartialPairCorrelationFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    partial_pair_correlation_function::T
    function PartialPairCorrelationFunction(radii::R, pcf::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, pcf)
        new{R,T,D,P}(radii, pcf)
    end
end

getargument(f::PartialPairCorrelationFunction) = f.radii
getestimate(f::PartialPairCorrelationFunction) = f.partial_pair_correlation_function
getextrafields(::PartialPairCorrelationFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function partial_paircorrelation_function(k::PartialKFunction{R,T,D}; penalty = 0.0, method = PCFMethodC()) where {R,T,D}
    pcf = Dict(index => K2paircorrelation(k.radii, val, Val{D}(), penalty, method) for (index, val) in k.partial_K_function)
    return PartialPairCorrelationFunction(k.radii, pcf, Val{D}())
end

function partial_paircorrelation_function(
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
    k = partial_K_function(
        data,
        region,
        radii,
        indices;
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
    )
    return partial_paircorrelation_function(k; penalty = penalty, method = method)
end

## direct method

function partial_paircorrelation_function_direct(
    f::PartialSpectra{D,F,P,N},
    λ,
    radii,
    indices = default_indices(f),
) where {D,F,P,N}
    check_spectral_indices(indices, f)
    pcf = Dict(
        index => sdf2pcf(f[index...], λ[index[1]], λ[index[2]], radii) for index in indices
    )
    return PartialPairCorrelationFunction(radii, pcf, Val{D}())
end

function partial_paircorrelation_function_direct(
    f::SpectralEstimate,
    λ,
    radii,
    indices::AbstractVector{Tuple{Int,Int}} = default_indices(f),
)
    partial_paircorrelation_function_direct(partial_spectra(f), λ, radii, indices)
end

function partial_paircorrelation_function_direct(
    f::SpectralEstimate{D,F,P,N},
    λ,
    radii,
    indices::AbstractVector{<:Tuple{Int,Int,AbstractVector{Int},AbstractVector{Int}}},
) where {D,F,P,N}
    check_spectral_indices(indices, f)
    C = Dict(
        index => sdf2pcf(
            partial_spectra(f, index[1], index[2], index[3], index[4]),
            λ[index[1]],
            λ[index[2]],
            radii,
        ) for index in indices
    )
    return PartialCFunction(radii, C, Val{D}())
end

function partial_paircorrelation_function_direct(
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
    return partial_paircorrelation_function_direct(f, λ, radii, indices)
end
