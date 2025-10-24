function apply_transform!(transform!, output::AbstractArray{<:SMatrix},
        input::AbstractArray{<:SMatrix}, ::MultipleTupleTrait, args...)
    for i in eachindex(output)
        output[i] = transform!(input[i], args...)
    end
    return output
end

function apply_transform!(transform!, output::AbstractArray{<:Number},
        input::AbstractArray{<:Number}, ::MultipleVectorTrait, args...)
    for i in CartesianIndices(size(output)[3:end])
        output[:, :, i] = transform!(input[:, :, i], args...)
    end
    return output
end

function apply_transform!(transform!, output::AbstractArray{<:Number},
        input::AbstractArray{<:Number}, ::SingleProcessTrait, args...)
    for i in eachindex(output)
        output[i] = transform!(input[i], args...)
    end
    return output
end
