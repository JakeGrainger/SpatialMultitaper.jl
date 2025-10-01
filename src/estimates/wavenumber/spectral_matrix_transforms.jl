function apply_transform!(
        transform, power::AbstractArray{<:SMatrix}, ::MultipleTupleTrait, args...)
    transform_wrapped(x) = transform(x, args...)
    map!(transform_wrapped, power, power)
    return power
end

function apply_transform!(
        transform, power::AbstractArray{<:Number}, ::MultipleVectorTrait, args...)
    for i in CartesianIndices(size(power)[3:end])
        power[:, :, i] = transform(power[:, :, i], args...)
    end
    return power
end

function apply_transform!(
        transform, power::AbstractArray{<:Number}, ::SingleProcessTrait, args...)
    transform_wrapped(x) = transform(x, args...)
    map!(transform_wrapped, power, power)
    return power
end
