function pair_correlation_function(data, region; kwargs...)
    pair_correlation_function(spatial_data(data, region); kwargs...)
end
function pair_correlation_function(data::SpatialData; pcf_method = PCFMethodC(), kwargs...)
    pair_correlation_function(k_function(data; kwargs...); pcf_method = pcf_method)
end
function pair_correlation_function(est::KFunction{E}; pcf_method = PCFMethodC()) where {E}
    radii = get_evaluation_points(est)
    value = k2paircorrelation(est, pcf_method)
    processinfo = get_process_information(est)
    estimationinfo = get_estimation_information(est)
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

# function pair_correlation_function_direct(f::Spectra{E}; radii) where {E}
#     value = sdf2pcf(f, radii)
#     processinfo = get_process_information(f)
#     estimationinfo = get_estimation_information(f)
#     return PairCorrelationFunction{E}(radii, value, processinfo, estimationinfo)
# end

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
