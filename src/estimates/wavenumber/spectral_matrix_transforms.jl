function apply_transform(
        transform, power::AbstractArray{<:SMatrix}, ::MultipleTupleTrait, args...)
    return transform.(power, Ref.(args)...)
end

function apply_transform(
        transform, power::AbstractArray{<:Number}, ::MultipleVectorTrait, args...)
    transform_wrapped(x) = transform(x, args...)
    mapslices(transform_wrapped, power, dims = (1, 2))
end

function apply_transform(
        transform, power::AbstractArray{<:Number}, ::SingleProcessTrait, args...)
    transform_wrapped(x) = transform(x, args...)
    map(transform_wrapped, power)
end
