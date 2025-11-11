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

function validate_radii(radii)
    @argcheck all(radii .>= 0)
    # side_length = sides(boundingbox(getregion(data)))
    # @argcheck all(radii .<= Meshes.ustrip(minimum(side_length)))
    nothing
end

process_radii(radii, ::SpatialData) = radii
function process_radii(::Nothing, data::SpatialData)
    region = getregion(data)
    short_side = Meshes.ustrip(minimum(sides(boundingbox(region))))
    return range(0, short_side / 3, length = 100)
end
function process_radii(::Nothing, ::AbstractEstimate)
    throw(ArgumentError("Default `radii` only available when computing from `SpatialData`, `radii` keyword argument must be provided."))
end

function process_pcf_radii(radii, data::SpatialData)
    return radii[(iszero(radii[1]) + 1):end]
end
function process_pcf_radii(::Nothing, arg)
    return process_radii(nothing, arg)[2:end]
end
function process_pcf_radii(radii, arg::KFunction)
    if iszero(radii[1])
        error("When computing pcf from a K function, you need to not use zero radii.")
    else
        return radii
    end
end
