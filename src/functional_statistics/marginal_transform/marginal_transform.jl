"""
    MarginallyTransformedEstimate{E, D, S, F, N, A, T, IP, IE} <: AbstractEstimate{E, D, N}

A wrapper for estimates that have undergone element-wise transformations.

This structure represents an estimate where a mathematical function has been applied
element-wise to the estimate values while preserving the argument structure and metadata.
Common examples include taking the real part, magnitude, phase, or logarithm of estimates.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `S`: Type of the original estimate before transformation
- `F`: Type of the transformation function
- `N`: Number of wavenumber dimensions
- `A`: Type of the argument (e.g., wavenumber grid)
- `T`: Type of the transformed estimate values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Fields
- `argument`: The argument (typically wavenumber) of the estimate
- `estimate`: The transformed estimate values
- `processinformation`: Information about the processes
- `estimationinformation`: Information about the estimation procedure

# Examples
```julia
# Take magnitude of a coherence estimate
mag_coh = abs(coherence_estimate)

# Extract real part of cross-spectrum
real_spec = real(cross_spectrum)

# Compute log-magnitude spectrum
log_spec = log(abs(spectrum))
```
"""
struct MarginallyTransformedEstimate{E, D, S, F, N, A, T, IP, IE} <:
       AbstractEstimate{E, D, N}
    argument::A
    estimate::T
    processinformation::IP
    estimationinformation::IE
    function MarginallyTransformedEstimate{E, S, N, F}(argument::A, estimate::T,
            processinfo::ProcessInformation{D}, estimationinfo::IE) where {
            E, S, N, F, A, T, D, IE}
        checkinputs(argument, estimate, processinfo)
        IP = typeof(processinfo)
        return new{E, D, S, F, N, A, T, IP, IE}(
            argument, estimate, processinfo, estimationinfo)
    end
end

## required interface

function computed_from(::Type{<:MarginallyTransformedEstimate{
        E, D, S, F}}) where {E, D, S, F}
    S
end

function allocate_estimate_memory(T::Type{<:MarginallyTransformedEstimate{E, D, S}},
        ::Type{<:S}, relevant_memory; kwargs...) where {E, D, S}
    transform = get_transform(T)
    possible_storetypes = Base.return_types(x -> transform.(x), (eltype(relevant_memory),))
    storetype = if length(possible_storetypes) == 1
        first(possible_storetypes)
    elseif length(possible_storetypes) > 1
        @warn "Transform is not type stable, so selecting Any as output type."
        Any
    else
        error("Transformation $transform has no valid return type for input of type $(eltype(relevant_memory)).")
    end
    out = zeros(storetype, size(relevant_memory))
    return out, nothing
end

function extract_relevant_memory(
        ::Type{<:MarginallyTransformedEstimate}, est::AbstractEstimate)
    get_estimates(est)
end
function extract_relevant_memory(
        ::Type{<:MarginallyTransformedEstimate}, mem::EstimateMemory)
    mem.output_memory
end

validate_core_parameters(::Type{<:MarginallyTransformedEstimate}; kwargs...) = nothing

function resolve_missing_parameters(::Type{<:MarginallyTransformedEstimate}, arg; kwargs...)
    kwargs
end

function validate_memory_compatibility(
        ::Type{T}, mem, arg; kwargs...) where {T <: MarginallyTransformedEstimate}
    @argcheck size(get_estimates(arg)) == size(mem.output_memory)
    transform = get_transform(T)
    @argcheck transform.(first(get_estimates(arg))) isa eltype(mem.output_memory)
    nothing
end

function compute_estimate!(
        ::Type{T}, mem, arg::AbstractEstimate{E, D, N};
        kwargs...) where {T <: MarginallyTransformedEstimate, E, D, N}
    transform = get_transform(T)
    estimate = apply_marginal_transform!(transform, mem.output_memory, get_estimates(arg))

    argument = get_evaluation_points(arg)
    processinfo = get_process_information(arg)
    estimationinfo = get_estimation_information(arg)
    S = typeof(arg)
    return MarginallyTransformedEstimate{E, S, N, typeof(transform)}(
        argument, estimate, processinfo, estimationinfo)
end

get_estimates(est::MarginallyTransformedEstimate) = est.estimate

get_evaluation_points(est::MarginallyTransformedEstimate) = est.argument

## additional interface

function get_estimate_name(T::Type{<:MarginallyTransformedEstimate})::String
    return "$(get_transform_name(T))($(get_estimate_name(getoriginaltype(T))))"
end
function get_short_estimate_name(T::Type{<:MarginallyTransformedEstimate})::String
    return "$(get_transform_name(T))($(get_short_estimate_name(getoriginaltype(T))))"
end
function get_short_base_estimate_name(T::Type{<:MarginallyTransformedEstimate})::String
    return "$(get_transform_name(T))($(get_short_base_estimate_name(getoriginaltype(T))))"
end

"""
    _construct_estimate_subset(::Type{<:MarginallyTransformedEstimate{...}}, trait, args...)

Internal constructor for creating estimate subsets with different traits.
"""
function _construct_estimate_subset(
        ::Type{<:MarginallyTransformedEstimate{E, D, S, F, N}},
        trait::Type{<:EstimateTrait},
        args...
)::MarginallyTransformedEstimate where {E, D, S, F, N}
    return MarginallyTransformedEstimate{trait, S, N, F}(args...)
end

# Type introspection utilities

"""
    getoriginaltype(::Type{<:MarginallyTransformedEstimate{E, D, N, S}}) where {E, D, N, S}

Extract the original estimate type before transformation.
"""
function getoriginaltype(::Type{<:MarginallyTransformedEstimate{
        E, D, S}}) where {E, D, S}
    return S
end

"""
    get_transform(est::MarginallyTransformedEstimate)

Get the transformation function applied to the estimate.
"""
get_transform(est::MarginallyTransformedEstimate) = get_transform(typeof(est))

function get_transform(::Type{<:MarginallyTransformedEstimate{
        E, D, S, F}}) where {E, D, S, F}
    return F.instance
end

"""
    get_transform_name(est::MarginallyTransformedEstimate)

Get a string representation of the transformation function name.
"""
get_transform_name(est::MarginallyTransformedEstimate) = get_transform_name(typeof(est))

function get_transform_name(T::Type{<:MarginallyTransformedEstimate})::String
    return string(get_transform(T))
end

## internals
"""
    apply_marginal_transform(transform::F, x::AbstractArray{<:SMatrix}) where {F}

Apply transformation to arrays of static matrices.

For arrays containing static matrices (common in multi-process spectral estimates),
applies the transformation element-wise to each matrix element.
"""
function apply_marginal_transform!(transform::F, mem, x::AbstractArray{<:SMatrix}) where {F}
    return map!(y -> transform.(y), mem, x)
end

"""
    apply_marginal_transform(transform, x::AbstractArray{<:Number})

Apply transformation to arrays of numbers.

For arrays of scalar values, applies the transformation element-wise across the array.
"""
function apply_marginal_transform!(transform::F, mem, x::AbstractArray{<:Number}) where {F}
    return map!(transform, mem, x)
end
