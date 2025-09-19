function apply_transform(
        transform::F,
        power::AbstractArray{SMatrix{P, P, T, L}, D},
        args...
) where {F, P, T, L, D}
    return transform.(power, Ref.(args)...)
end

function apply_transform(transform, power::AbstractArray{<:Number, D}, args...) where {D}
    mapslices(transform, power, dims = (1, 2), args...)
end
