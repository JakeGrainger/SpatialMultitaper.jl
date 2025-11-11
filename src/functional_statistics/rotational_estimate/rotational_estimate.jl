"""
    RotationalEstimate{E, D, S, A, T, IP, IE} <: IsotropicEstimate{E, D}

An estimate that has undergone rotational averaging to produce isotropic results.

Rotational averaging transforms anisotropic estimates (which vary with direction)
into isotropic estimates that depend only on radial distance from the origin.
This is commonly used in spectral analysis to study scale-dependent properties
independent of orientation.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `S`: Type of the original anisotropic estimate
- `A`: Type of radii array
- `T`: Type of the rotationally averaged estimate values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Fields
- `radii::A`: Radial distances at which the estimate is evaluated
- `estimate::T`: Rotationally averaged estimate values
- `processinformation::IP`: Information about the processes
- `estimationinformation::IE`: Information about the estimation procedure

# Examples
```julia
# Create rotational estimate with default parameters
rot_spec = rotational_estimate(anisotropic_spectrum)

# Create with custom radii and kernel
radii = 0.1:0.1:1.0
kernel = GaussKernel(0.05)
rot_spec = rotational_estimate(spectrum, radii=radii, kernel=kernel)
```
"""
struct RotationalEstimate{E, D, S, A, T, IP, IE} <: IsotropicEstimate{E, D}
    radii::A
    estimate::T
    processinformation::IP
    estimationinformation::IE
    function RotationalEstimate{E, S}(
            radii::A, estimate::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, S, A, T, D, IE}
        checkinputs(radii, estimate, processinfo)
        IP = typeof(processinfo)
        return new{E, D, S, A, T, IP, IE}(
            radii, estimate, processinfo, estimationinfo)
    end
end

## interface

computed_from(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S} = S

function allocate_estimate_memory(::Type{<:RotationalEstimate{E, D, S}},
        ::Type{<:S}, relevant_memory; kwargs...) where {E, D, S}
    out = preallocate_rotational_radial_output(relevant_memory...; kwargs...)
    return out, nothing
end

function preallocate_rotational_radial_output(
        source, ::Union{SingleProcessTrait, MultipleTupleTrait}; rotation_radii, kwargs...)
    return zeros(real(eltype(source)), length(rotation_radii))
end
function preallocate_rotational_radial_output(
        source, ::MultipleVectorTrait; rotation_radii, kwargs...)
    return zeros(
        real(eltype(source)), size(source, 1), size(source, 2), length(rotation_radii))
end

function extract_relevant_memory(::Type{<:RotationalEstimate}, source::AbstractEstimate)
    return get_estimates(source), process_trait(source)
end

function extract_relevant_memory(::Type{<:RotationalEstimate}, source::EstimateMemory)
    return source.output_memory, process_trait(source)
end

function validate_core_parameters(::Type{<:RotationalEstimate}; kwargs...)
    if haskey(kwargs, :rotation_radii)
        rotation_radii = kwargs[:rotation_radii]
        @argcheck rotation_radii isa AbstractVector{<:Real}
        @argcheck all(r -> r ≥ 0, rotation_radii)
    end
    if haskey(kwargs, :kernel)
        kernel = kwargs[:kernel]
        @argcheck kernel isa RotationalKernel
    end
end

function resolve_missing_parameters(::Type{<:RotationalEstimate}, arg; kwargs...)
    rotation_radii = if haskey(kwargs, :rotation_radii)
        kwargs[:rotation_radii]
    else
        default_rotational_radii(arg; kwargs...)
    end
    kernel = if haskey(kwargs, :kernel)
        kwargs[:kernel]
    else
        default_rotational_kernel(arg; rotation_radii = rotation_radii, kwargs...)
    end
    other_kwargs = filter(kv -> kv[1] ∉ (:rotation_radii, :kernel), kwargs)
    return (; rotation_radii = rotation_radii, kernel = kernel, other_kwargs...)
end

function validate_memory_compatibility(
        ::Type{<:RotationalEstimate}, mem, source; rotation_radii, kwargs...)
    validate_radial_memory(mem.output_memory, process_trait(source), rotation_radii)
end

function compute_estimate!(::Type{<:RotationalEstimate{E, D, S}}, mem,
        source::S; rotation_radii, kernel, kwargs...) where {E, D, S}
    x = get_evaluation_points(source)
    y = get_estimates(source)
    _smoothed_rotational!(
        mem.output_memory, x, y, process_trait(source), rotation_radii, kernel)

    processinfo = get_process_information(source)
    estimationinfo = get_estimation_information(source)
    return RotationalEstimate{E, S}(
        rotation_radii, mem.output_memory, processinfo, estimationinfo)
end

get_evaluation_points(f::RotationalEstimate) = f.radii

get_estimates(f::RotationalEstimate) = f.estimate

## additional interface

"""
    getoriginaltype(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S}

Extract the original estimate type before rotational averaging.
"""
getoriginaltype(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S} = S

"""
    get_estimate_name(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S}

Get the name of the estimate, reflecting its rotational and original traits.

The name indicates the order of application of rotational and partial traits.
"""
function get_estimate_name(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S}
    original_name = get_base_estimate_name(S)

    if E === MarginalTrait
        return "rotational " * original_name
    elseif E === PartialTrait
        # Check if original was marginal or partial to determine order
        original_trait = _extract_estimate_trait(S)
        if original_trait === MarginalTrait
            return "partial rotational " * original_name  # Partial applied after rotational
        else
            return "rotational partial " * original_name  # Rotational applied after partial
        end
    end
end
function get_short_estimate_name(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S}
    original_name = get_short_base_estimate_name(S)

    if E === MarginalTrait
        return "rotational " * original_name
    elseif E === PartialTrait
        # Check if original was marginal or partial to determine order
        original_trait = _extract_estimate_trait(S)
        if original_trait === MarginalTrait
            return "partial rotational " * original_name  # Partial applied after rotational
        else
            return "rotational partial " * original_name  # Rotational applied after partial
        end
    end
end
function get_short_base_estimate_name(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S}
    original_name = get_short_base_estimate_name(S)

    if E === MarginalTrait
        return "rotational " * original_name
    elseif E === PartialTrait
        # Check if original was marginal or partial to determine order
        original_trait = _extract_estimate_trait(S)
        if original_trait === MarginalTrait
            return "partial rotational " * original_name  # Partial applied after rotational
        else
            return "rotational partial " * original_name  # Rotational applied after partial
        end
    end
end

# Helper function to extract trait from type
_extract_estimate_trait(::Type{<:AbstractEstimate{T}}) where {T} = T

"""
    _construct_estimate_subset(
        ::Type{<:RotationalEstimate{E, D, S}},
        trait::Type{<:EstimateTrait},
        argument, estimate, processinfo, estimationinfo
    ) where {E, D, S}

Construct a subset estimate with a different trait from an existing RotationalEstimate.

This is used to create partial or marginal estimates from a rotational estimate.
"""
function _construct_estimate_subset(
        ::Type{<:RotationalEstimate{E, D, S}},
        trait::Type{<:EstimateTrait},
        argument, estimate, processinfo, estimationinfo
) where {E, D, S}
    # For RotationalEstimate, we need {E, S} constructor
    return RotationalEstimate{trait, S}(
        argument, estimate, processinfo, estimationinfo
    )
end

## internals
include("rotational_kernels.jl")

function _smoothed_rotational!(out, x::NTuple{D}, y::AbstractArray{T, D},
        trait::Union{SingleProcessTrait, MultipleTupleTrait}, rotation_radii,
        kernel) where {D, T <: Union{<:Number, <:SMatrix}}
    @argcheck length(x) == ndims(y)
    @argcheck size(y) == length.(x)
    @argcheck length(out) == length(rotation_radii)
    xitr = Iterators.ProductIterator(x)
    for (i, r) in enumerate(rotation_radii)
        num = sum(f * kernel(norm(u) - r) for (u, f) in zip(xitr, y))
        denom = sum(kernel(norm(u) - r) for u in xitr)
        out[i] = real(num / denom)
    end
    return out
end

function _smoothed_rotational!(
        out::AbstractArray{<:Number, 3}, x::NTuple{D}, y::AbstractArray{<:Number, N},
        ::MultipleVectorTrait, rotation_radii, kernel) where {D, N}
    @argcheck length(x) <= ndims(y)
    @argcheck size(y)[3:end] == length.(x)

    _trait = SingleProcessTrait()
    @views for idx in CartesianIndices(size(y)[1:2])
        y_slice = y[idx, ntuple(Returns(:), Val{N - 2}())...]
        _smoothed_rotational!(out[idx, :], x, y_slice, _trait, rotation_radii, kernel)
    end
    return out
end

## defaults
function default_rotational_radii(::SpatialData; nk, kmax, kwargs...)
    return default_rotational_radii(_choose_wavenumbers_1d.(nk, kmax))
end

function default_rotational_radii(est::AbstractEstimate; kwargs...)
    return default_rotational_radii(get_evaluation_points(est))
end

function default_rotational_radii(wavenumber::NTuple{D, AbstractVector{<:Real}}) where {D}
    max_wavenumber = minimum(x -> maximum(abs, x), wavenumber)
    step_wavenumber = minimum(step, wavenumber)
    # Offset by half step to integrate with endpoints at zero range steps
    used_range = range(step_wavenumber / 2, stop = max_wavenumber, step = step_wavenumber)
    return used_range
end

function default_rotational_kernel(::SpatialData; nk, kmax, kwargs...)
    return default_rotational_kernel(_choose_wavenumbers_1d.(nk, kmax))
end

function default_rotational_kernel(est::AbstractEstimate; kwargs...)
    return default_rotational_kernel(get_evaluation_points(est))
end

function default_rotational_kernel(x::NTuple{D, AbstractVector{<:Real}}) where {D}
    max_step = maximum(step, x) # TODO: currently only supports uniform steps
    return RectKernel(2 * max_step)
end

function default_rotational_kernel(x::AbstractVector{<:Real})
    max_step = step(x)
    return RectKernel(2 * max_step)
end
