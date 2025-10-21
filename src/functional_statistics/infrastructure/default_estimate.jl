function compute(::Type{T}, arg; kwargs...) where {T <: AbstractEstimate}
    resolved_kwargs = resolve_parameters(T, arg; kwargs...)
    mem = preallocate_memory(T, arg; resolved_kwargs...)
    return estimate_function!(T, mem, arg; resolved_kwargs...)
end

function resolve_parameters(::Type{T}, arg; kwargs...) where {T <: AbstractEstimate}
    validate_core_parameters(T; kwargs...)
    # Resolve parameters that may depend on the input argument's properties
    resolved = resolve_missing_parameters(T; kwargs...)
    return apply_parameter_defaults(T; resolved...)
end

function estimate_function!(
        ::Type{T}, mem::EstimateMemory{T}, arg; kwargs...) where {T <: AbstractEstimate}
    check_memory_compatibility(T, mem, arg; kwargs...)
    estimate_function!(computed_from(T), mem.previous_memory, arg; kwargs...)
    return compute_estimate!(T, mem, mem.previous_memory; kwargs...)
end

function estimate_function!(
        ::Type{T}, mem::BaseEstimateMemory{T}, arg; kwargs...) where {T <: AbstractEstimate}
    check_memory_compatibility(T, mem, arg; kwargs...)
    return compute_estimate!(T, mem, arg; kwargs...)
end

struct EstimateMemory{T, M, S <: Union{Nothing, <:EstimateMemory}}
    estimate_type::Type{T}
    memory_data::M
    previous_memory::S
end
const BaseEstimateMemory{T, M} = EstimateMemory{T, M, Nothing}

function preallocate_memory(::Type{T}, arg; kwargs...) where {T <: AbstractEstimate}
    validate_computation_chain(T, arg)
    return _preallocate_memory(T, computed_from(T), arg; kwargs...)
end

function _preallocate_memory(::Type{T}, ::Type{S}, arg; kwargs...) where {T, S}
    previous_memory = _preallocate_memory(S, computed_from(S), arg; kwargs...)
    memory_data = allocate_estimate_memory(T, S, previous_memory; kwargs...)
    return EstimateMemory(T, memory_data, previous_memory)
end
function _preallocate_memory(::Type{T}, ::Type{S}, arg::S; kwargs...) where {T, S}
    # Extract the memory structure that was used to compute this existing estimate
    previous_memory = extract_allocation_memory(arg)
    memory_data = allocate_estimate_memory(T, S, previous_memory; kwargs...)
    return EstimateMemory(T, memory_data, previous_memory)
end

function validate_computation_chain(T, arg)
    S = computed_from(T)
    if arg isa S
        return nothing
    elseif S === SpatialData # If you reach this point without typeof(arg) in the chain, you cant go further as Spatial Data is the base of everything
        error("Estimates of type $T cannot be computed from $(typeof(arg)).")
    elseif S === T
        error("Estimates of type $T have an incorrectly defined computation chain, you will have missdefined `computed_from` to be circular.")
    else
        validate_computation_chain(S, arg)
    end
end
