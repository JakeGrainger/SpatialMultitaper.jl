## spectra

function spectra(data, region::Meshes.Geometry; kwargs...)
    return spectra(spatial_data(data, region); kwargs...)
end

function spectra(data::SpatialData; kwargs...)
    return compute(Spectra{MarginalTrait}, data; kwargs...)
end

function partial_spectra(data, region::Meshes.Geometry; kwargs...)
    return partial_spectra(spatial_data(data, region); kwargs...)
end

function partial_spectra(arg; kwargs...)
    return compute(Spectra{PartialSpectra}, arg; kwargs...)
end

## coherence

function coherence(data, region::Meshes.Geometry; kwargs...)
    return coherence(spatial_data(data, region); kwargs...)
end

function coherence(arg; kwargs...)
    return compute(Coherence{MarginalTrait}, arg; kwargs...)
end

function coherence(arg::PartialAbstractEstimate; kwargs...)
    return compute(Coherence{PartialSpectra}, arg; kwargs...)
end

function partial_coherence(data, region::Meshes.Geometry; kwargs...)
    return partial_coherence(spatial_data(data, region); kwargs...)
end

function partial_coherence(arg; kwargs...)
    return compute(Coherence{PartialSpectra}, arg; kwargs...)
end
