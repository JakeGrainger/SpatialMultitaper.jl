function preallocate_spatial_output(spectrum::NormalOrRotationalSpectra, radii)
    return _preallocate_spatial_output(
        get_estimates(spectrum), process_trait(spectrum), radii)
end
function _preallocate_spatial_output(
        power, ::Union{SingleProcessTrait, MultipleTupleTrait}, radii)
    return zeros(real(eltype(power)), length(radii))
end
function _preallocate_spatial_output(power, ::MultipleVectorTrait, radii)
    return zeros(real(eltype(power)), size(power, 1), size(power, 2), length(radii))
end
