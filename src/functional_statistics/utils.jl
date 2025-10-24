function preallocate_radial_output(
        source, ::Union{SingleProcessTrait, MultipleTupleTrait}; radii, kwargs...)
    return zeros(real(eltype(source)), length(radii))
end
function preallocate_radial_output(source, ::MultipleVectorTrait; radii, kwargs...)
    return zeros(real(eltype(source)), size(source, 1), size(source, 2), length(radii))
end

function validate_radial_memory(
        mem::AbstractVector, ::Union{SingleProcessTrait, MultipleTupleTrait}, radii)
    @argcheck length(mem) == length(radii)
end
function validate_radial_memory(mem::AbstractArray, ::MultipleVectorTrait, radii)
    @argcheck size(mem, 3) == length(radii)
    @argcheck ndims(mem) == 3
end
