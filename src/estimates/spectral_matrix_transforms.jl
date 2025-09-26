function apply_transform(
        transform::F,
        freq,
        power::AbstractArray{<:SMatrix, D},
        args...
) where {F, D}
    return transform.(power, Ref.(args)...)
end

function apply_transform(
        transform, freq, power::AbstractArray{<:Number, D}, args...) where {D}
    transform_wrapped(x) = transform(x, args...)
    mapslices(transform_wrapped, power, dims = (1, 2))
end

function apply_transform(
        transform, freq::NTuple{D}, power::AbstractArray{<:Number, D}, args...) where {D}
    transform_wrapped(x) = transform(x, args...)
    map(transform_wrapped, power)
end
