function apply_transform(
        transform::F,
        power::AbstractArray{SMatrix{P, P, T, L}, D},
        args...
) where {F, P, T, L, D}
    return transform.(power, Ref.(args)...)
end

function apply_transform(transform, power::AbstractArray{<:Number, D}, args...) where {D}
    transform_wrapped(x) = transform(x, args...)
    mapslices(transform_wrapped, power, dims = (1, 2))
end
