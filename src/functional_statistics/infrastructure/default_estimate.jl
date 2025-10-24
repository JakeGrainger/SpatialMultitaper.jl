struct EstimateMemory{T, M, I, S <: Union{Nothing, <:EstimateMemory}}
    estimate_type::Type{T}
    output_memory::M
    internal_memory::I
    previous_memory::S
end
const BaseEstimateMemory{T, M, I} = EstimateMemory{T, M, I, Nothing}

function compute(::Type{T}, arg; kwargs...) where {T <: AbstractEstimate}
    kwargs = filter(kv -> kv[2] !== nothing, kwargs) # required for the api and R implementation
    resolved_kwargs = resolve_parameters(T, arg; kwargs...)
    mem = preallocate_memory(T, arg; resolved_kwargs...)
    return estimate_function!(T, mem, arg; resolved_kwargs...)
end

function resolve_parameters(::Type{T}, arg; kwargs...)
    previous_resolved = resolve_parameters(computed_from(T), arg; kwargs...)

    # Phase 1: Validate what the user actually provided
    validate_core_parameters(T; previous_resolved...)

    # Phase 2: Let each estimate type handle its own parameter resolution + defaults
    # This is where Spectra handles dk/nk/kmax constraints, CFunction handles radii, etc.
    resolved = resolve_missing_parameters(T, arg; previous_resolved...)

    return resolved
end

function resolve_parameters(::Type{T}, arg::T; kwargs...) where {T}
    return kwargs # base case
end

function estimate_function!(
        ::Type{T}, mem::EstimateMemory{T}, arg; kwargs...) where {T <: AbstractEstimate}
    check_memory_compatibility(T, mem, arg; kwargs...)
    previous_arg = estimate_function!(computed_from(T), mem.previous_memory, arg; kwargs...)
    return compute_estimate!(T, mem, previous_arg; kwargs...)
end

function estimate_function!(
        ::Type{T}, mem::BaseEstimateMemory{T}, arg; kwargs...) where {T <: AbstractEstimate}
    check_memory_compatibility(T, mem, arg; kwargs...)
    return compute_estimate!(T, mem, arg; kwargs...)
end

function check_memory_compatibility(
        ::Type{T}, mem::EstimateMemory{T}, arg; kwargs...) where {T <: AbstractEstimate}
    # Memory type compatibility is guaranteed by type system
    # Type-specific compatibility checks
    validate_memory_compatibility(T, mem, arg; kwargs...)
    return nothing
end

function preallocate_memory(::Type{T}, arg; kwargs...) where {T <: AbstractEstimate}
    validate_computation_chain(T, arg)
    source_type = select_source_type(T, arg)
    return _preallocate_memory(T, source_type, arg; kwargs...)
end

function select_source_type(::Type{T}, arg) where {T}
    possible_sources = computed_from(T)

    # Handle single type
    if possible_sources isa Type
        return possible_sources
    end

    # Handle tuple - select the one that matches arg type
    for S in possible_sources
        if arg isa S
            return S
        end
    end

    # If arg doesn't match any source directly, use the preferred (first) source
    # This handles cases where we need to build a computation chain
    return first(possible_sources)
end

function _preallocate_memory(::Type{T}, ::Type{S}, arg; kwargs...) where {T, S}
    previous_memory = _preallocate_memory(S, select_source_type(S, arg), arg; kwargs...)
    relevant_memory = extract_relevant_memory(T, previous_memory)
    output_memory, internal_memory = allocate_estimate_memory(
        T, S, relevant_memory; kwargs...)
    return EstimateMemory(T, output_memory, internal_memory, previous_memory)
end
function _preallocate_memory(::Type{T}, ::Type{S}, arg::S; kwargs...) where {T, S}
    # Extract the memory structure that was used to compute this existing estimate
    relevant_memory = extract_relevant_memory(T, arg) # only needed to decide the required memory
    output_memory, internal_memory = allocate_estimate_memory(
        T, S, relevant_memory; kwargs...)
    return EstimateMemory(T, output_memory, internal_memory, nothing)
end

function validate_computation_chain(T, arg)
    possible_sources = computed_from(T)
    possible_sources = possible_sources isa Type ? (possible_sources,) : possible_sources
    for S in possible_sources
        if arg isa S
            return nothing  # Valid path found
        end
    end
    if S === SpatialData # If you reach this point without typeof(arg) in the chain, you cant go further as Spatial Data is the base of everything
        error("Estimates of type $T cannot be computed from $(typeof(arg)).")
    elseif S === T
        # this error will only happen if an incorrect extension is made when adding a new statistic
        error("Estimates of type $T have an incorrectly defined computation chain, you will have missdefined `computed_from` to be circular.")
    else
        validate_computation_chain(S, arg)
    end
end
