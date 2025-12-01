## spectra

function spectra(data, region::Meshes.Geometry; kwargs...)
    return spectra(spatial_data(data, region); kwargs...)
end

functional_statistic_type(::typeof(spectra), arg) = Spectra{MarginalTrait, embeddim(arg)}

function spectra(arg; kwargs...)
    return compute(functional_statistic_type(spectra, arg), arg; kwargs...)
end

function partial_spectra(data, region::Meshes.Geometry; kwargs...)
    return partial_spectra(spatial_data(data, region); kwargs...)
end

function functional_statistic_type(::typeof(partial_spectra), arg)
    Spectra{PartialTrait, embeddim(arg)}
end

function partial_spectra(arg; kwargs...)
    return compute(functional_statistic_type(partial_spectra, arg), arg; kwargs...)
end

## coherence

function coherence(data, region::Meshes.Geometry; kwargs...)
    return coherence(spatial_data(data, region); kwargs...)
end

function functional_statistic_type(::typeof(coherence), arg)
    Coherence{MarginalTrait, embeddim(arg)}
end

function coherence(arg; kwargs...)
    return compute(functional_statistic_type(coherence, arg), arg; kwargs...)
end

function coherence(arg::PartialAbstractEstimate; kwargs...)
    return compute(Coherence{PartialTrait, embeddim(arg)}, arg; kwargs...)
end

function partial_coherence(data, region::Meshes.Geometry; kwargs...)
    return partial_coherence(spatial_data(data, region); kwargs...)
end

function functional_statistic_type(::typeof(partial_coherence), arg)
    Coherence{PartialTrait, embeddim(arg)}
end

function partial_coherence(arg; kwargs...)
    return compute(functional_statistic_type(partial_coherence, arg), arg; kwargs...)
end

## uncorrected partial spectra

function partial_spectra_uncorrected(arg::SpatialData; kwargs...)
    spectrum = spectra(arg; kwargs...)
    # Create a modified spectrum with no taper information for uncorrected computation
    new_spectrum = Spectra{MarginalTrait}(
        get_evaluation_points(spectrum), get_estimates(spectrum),
        get_process_information(spectrum), EstimationInformation(nothing))
    return partial_spectra(new_spectrum)
end

## other transforms

magnitude_coherence(arg; kwargs...) = abs(coherence(arg; kwargs...))
function magnitude_coherence(data, region; kwargs...)
    magnitude_coherence(spatial_data(data, region); kwargs...)
end

function partial_magnitude_coherence(arg; kwargs...)
    abs(partial_coherence(arg; kwargs...))
end
function partial_magnitude_coherence(data, region; kwargs...)
    partial_magnitude_coherence(spatial_data(data, region); kwargs...)
end

magnitude_squared_coherence(arg; kwargs...) = abs2(coherence(arg; kwargs...))
function magnitude_squared_coherence(data, region; kwargs...)
    magnitude_squared_coherence(spatial_data(data, region); kwargs...)
end

function partial_magnitude_squared_coherence(arg; kwargs...)
    return abs2(partial_coherence(arg; kwargs...))
end
function partial_magnitude_squared_coherence(data, region; kwargs...)
    return partial_magnitude_squared_coherence(spatial_data(data, region); kwargs...)
end

phase(arg; kwargs...) = angle(coherence(arg; kwargs...))
phase(data, region; kwargs...) = phase(spatial_data(data, region); kwargs...)

partial_phase(arg; kwargs...) = angle(partial_coherence(arg; kwargs...))
function partial_phase(data, region; kwargs...)
    partial_phase(spatial_data(data, region); kwargs...)
end
