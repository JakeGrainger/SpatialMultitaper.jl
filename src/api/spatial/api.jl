include("docstrings.jl")

## CFunction
function functional_statistic_type(::typeof(c_function), arg)
    CFunction{MarginalTrait, embeddim(arg)}
end

function c_function(data, region; kwargs...)
    return c_function(spatial_data(data, region); kwargs...)
end

function c_function(arg; kwargs...)
    compute(functional_statistic_type(c_function, arg), arg; kwargs...)
end

function c_function(arg::PartialAbstractEstimate; kwargs...)
    partial_c_function(arg; kwargs...)
end

function functional_statistic_type(::typeof(partial_c_function), arg)
    CFunction{PartialTrait, embeddim(arg)}
end

function partial_c_function(data, region; kwargs...)
    return partial_c_function(spatial_data(data, region); kwargs...)
end

function partial_c_function(arg; kwargs...)
    compute(functional_statistic_type(partial_c_function, arg), arg; kwargs...)
end

## KFunction
function functional_statistic_type(::typeof(k_function), arg)
    KFunction{MarginalTrait, embeddim(arg)}
end

function k_function(data, region; kwargs...)
    return k_function(spatial_data(data, region); kwargs...)
end

function k_function(arg; kwargs...)
    compute(functional_statistic_type(k_function, arg), arg; kwargs...)
end

function k_function(arg::PartialAbstractEstimate; kwargs...)
    partial_k_function(arg; kwargs...)
end

function functional_statistic_type(::typeof(partial_k_function), arg)
    KFunction{PartialTrait, embeddim(arg)}
end

function partial_k_function(data, region; kwargs...)
    return partial_k_function(spatial_data(data, region); kwargs...)
end

function partial_k_function(arg; kwargs...)
    compute(functional_statistic_type(partial_k_function, arg), arg; kwargs...)
end

## LFunction
function functional_statistic_type(::typeof(l_function), arg)
    LFunction{MarginalTrait, embeddim(arg)}
end

function l_function(data, region; kwargs...)
    return l_function(spatial_data(data, region); kwargs...)
end

function l_function(arg; kwargs...)
    compute(functional_statistic_type(l_function, arg), arg; kwargs...)
end

function l_function(arg::PartialAbstractEstimate; kwargs...)
    partial_l_function(arg; kwargs...)
end

function functional_statistic_type(::typeof(partial_l_function), arg)
    LFunction{PartialTrait, embeddim(arg)}
end

function partial_l_function(data, region; kwargs...)
    return partial_l_function(spatial_data(data, region); kwargs...)
end

function partial_l_function(arg; kwargs...)
    compute(functional_statistic_type(partial_l_function, arg), arg; kwargs...)
end

## CenteredLFunction
function functional_statistic_type(::typeof(centered_l_function), arg)
    CenteredLFunction{MarginalTrait, embeddim(arg)}
end

function centered_l_function(data, region; kwargs...)
    return centered_l_function(spatial_data(data, region); kwargs...)
end

function centered_l_function(arg; kwargs...)
    compute(functional_statistic_type(centered_l_function, arg), arg; kwargs...)
end

function centered_l_function(arg::PartialAbstractEstimate; kwargs...)
    partial_centered_l_function(arg; kwargs...)
end

function functional_statistic_type(::typeof(partial_centered_l_function), arg)
    CenteredLFunction{PartialTrait, embeddim(arg)}
end

function partial_centered_l_function(data, region; kwargs...)
    return partial_centered_l_function(spatial_data(data, region); kwargs...)
end

function partial_centered_l_function(arg; kwargs...)
    compute(functional_statistic_type(partial_centered_l_function, arg), arg; kwargs...)
end

## Pair correlation function
function functional_statistic_type(::typeof(pair_correlation_function), arg)
    PairCorrelationFunction{MarginalTrait, embeddim(arg)}
end

function pair_correlation_function(data, region; kwargs...)
    return pair_correlation_function(spatial_data(data, region); kwargs...)
end
function pair_correlation_function(arg; kwargs...)
    compute(functional_statistic_type(pair_correlation_function, arg), arg; kwargs...)
end
function pair_correlation_function(arg::PartialAbstractEstimate; kwargs...)
    partial_pair_correlation_function(arg; kwargs...)
end
function functional_statistic_type(::typeof(partial_pair_correlation_function), arg)
    PairCorrelationFunction{PartialTrait, embeddim(arg)}
end

function partial_pair_correlation_function(data, region; kwargs...)
    return partial_pair_correlation_function(spatial_data(data, region); kwargs...)
end
function partial_pair_correlation_function(arg; kwargs...)
    compute(
        functional_statistic_type(partial_pair_correlation_function, arg), arg; kwargs...)
end
