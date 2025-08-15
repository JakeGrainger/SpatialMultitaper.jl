struct PartialPairCorrelationFunction{R,P}
    radii::R
    pcf::P
end

function partial_paircorrelation_function(
    f::PartialSpectra,
    λ,
    radii,
    indices = default_indices(f),
)
    check_spectral_indices(indices, f)
    pcf = Dict(
        index => sdf2pcf(f[index...], λ[index[1]], λ[index[2]], radii) for index in indices
    )
    return PartialPairCorrelationFunction(radii, pcf)
end

function partial_paircorrelation_function(
    f::SpectralEstimate,
    λ,
    radii,
    indices::AbstractVector{Tuple{Int,Int}} = default_indices(f),
)
    partial_paircorrelation_function(partial_spectra(f), λ, radii, indices)
end

function partial_paircorrelation_function(
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

function partial_paircorrelation_function(
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
    return partial_paircorrelation_function(f, λ, radii, indices)
end
