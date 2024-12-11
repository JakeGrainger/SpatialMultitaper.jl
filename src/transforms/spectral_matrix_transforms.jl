function apply_transform(transform::F, power::AbstractArray{SMatrix{P, P, T, L}, D}) where {F, P, T, L, D}
    return transform.(power)
end
function apply_transform(transform, power)
    mapslices(transform, power, dims = (1, 2))
end
