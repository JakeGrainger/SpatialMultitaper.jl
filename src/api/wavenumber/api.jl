## spectra

function spectra(data, region::Meshes.Geometry; kwargs...)
    return spectra(spatial_data(data, region); kwargs...)
end

function spectra(data::SpatialData; kwargs...)
    return compute(Spectra{MarginalTrait, embeddim(data)}, data; kwargs...)
end

function partial_spectra(data, region::Meshes.Geometry; kwargs...)
    return partial_spectra(spatial_data(data, region); kwargs...)
end

function partial_spectra(arg; kwargs...)
    return compute(Spectra{PartialSpectra, embeddim(data)}, arg; kwargs...)
end

## coherence

function coherence(data, region::Meshes.Geometry; kwargs...)
    return coherence(spatial_data(data, region); kwargs...)
end

function coherence(arg; kwargs...)
    return compute(Coherence{MarginalTrait, embeddim(arg)}, arg; kwargs...)
end

function coherence(arg::PartialAbstractEstimate; kwargs...)
    return compute(Coherence{PartialTrait, embeddim(arg)}, arg; kwargs...)
end

function partial_coherence(data, region::Meshes.Geometry; kwargs...)
    return partial_coherence(spatial_data(data, region); kwargs...)
end

function partial_coherence(arg; kwargs...)
    return compute(Coherence{PartialTrait, embeddim(arg)}, arg; kwargs...)
end

## uncorrected partial spectra

function partial_spectra_uncorrected(arg; kwargs...)
    spectrum = spectra(arg; kwargs...)
    # Create a modified spectrum with no taper information for uncorrected computation
    new_spectrum = Spectra{MarginalTrait}(
        get_evaluation_points(spectrum), get_estimates(spectrum),
        get_process_information(spectrum), EstimationInformation(nothing))
    return partial_spectra(new_spectrum)
end
