function apply_transform!(
        transform!, power::AbstractArray{<:SMatrix}, ::MultipleTupleTrait, args...)
    for i in eachindex(power)
        power[i] = transform!(power[i], args...)
    end
    return power
end

function apply_transform!(
        transform!, power::AbstractArray{<:Number}, ::MultipleVectorTrait, args...)
    for i in CartesianIndices(size(power)[3:end])
        power[:, :, i] = transform!(power[:, :, i], args...)
    end
    return power
end

function apply_transform!(
        transform!, power::AbstractArray{<:Number}, ::SingleProcessTrait, args...)
    for i in eachindex(power)
        power[i] = transform!(power[i], args...)
    end
    return power
end
