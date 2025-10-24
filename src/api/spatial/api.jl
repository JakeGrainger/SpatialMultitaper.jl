include("docstrings.jl")

## CFunction

function c_function(data, region; kwargs...)
    return c_function(spatial_data(data, region); kwargs...)
end

function c_function(arg; kwargs...)
    compute(CFunction{MarginalTrait}, arg; kwargs...)
end

function c_function(arg::PartialAbstractEstimate; kwargs...)
    partial_c_function(arg; kwargs...)
end

function partial_c_function(data, region; kwargs...)
    return partial_c_function(spatial_data(data, region); kwargs...)
end

function partial_c_function(arg; kwargs...)
    compute(CFunction{PartialTrait}, arg; kwargs...)
end

## KFunction

function k_function(data, region; kwargs...)
    return k_function(spatial_data(data, region); kwargs...)
end

function k_function(arg; kwargs...)
    compute(KFunction{MarginalTrait}, arg; kwargs...)
end

function k_function(arg::PartialAbstractEstimate; kwargs...)
    partial_k_function(arg; kwargs...)
end

function partial_k_function(data, region; kwargs...)
    return partial_k_function(spatial_data(data, region); kwargs...)
end

function partial_k_function(arg; kwargs...)
    compute(KFunction{PartialTrait}, arg; kwargs...)
end

## LFunction

function l_function(data, region; kwargs...)
    return l_function(spatial_data(data, region); kwargs...)
end

function l_function(arg; kwargs...)
    compute(LFunction{MarginalTrait}, arg; kwargs...)
end

function l_function(arg::PartialAbstractEstimate; kwargs...)
    partial_l_function(arg; kwargs...)
end

function partial_l_function(data, region; kwargs...)
    return partial_l_function(spatial_data(data, region); kwargs...)
end

function partial_l_function(arg; kwargs...)
    compute(LFunction{PartialTrait}, arg; kwargs...)
end

## CenteredLFunction

function centered_l_function(data, region; kwargs...)
    return centered_l_function(spatial_data(data, region); kwargs...)
end

function centered_l_function(arg; kwargs...)
    compute(CenteredLFunction{MarginalTrait}, arg; kwargs...)
end

function centered_l_function(arg::PartialAbstractEstimate; kwargs...)
    partial_centered_l_function(arg; kwargs...)
end

function partial_centered_l_function(data, region; kwargs...)
    return partial_centered_l_function(spatial_data(data, region); kwargs...)
end

function partial_centered_l_function(arg; kwargs...)
    compute(CenteredLFunction{PartialTrait}, arg; kwargs...)
end
