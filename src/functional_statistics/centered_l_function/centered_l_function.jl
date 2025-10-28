"""
    CenteredLFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}

Centered L function estimate for spatial pattern analysis.

The centered L function is defined as ``CL(r) = L(r) - r``, which removes the
theoretical expectation for a Poisson process. This makes it easier to detect
deviations from complete spatial randomness:

- ``CL(r) = 0`` indicates Poisson behavior (complete spatial randomness)
- ``CL(r) > 0`` indicates clustering at scale r
- ``CL(r) < 0`` indicates regularity/inhibition at scale r

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `A`: Type of radii array
- `T`: Type of centered L function values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Mathematical Background
For a Poisson process, ``L(r) = r`` theoretically. The centering transformation:
``CL(r) = L(r) - r`` removes this trend, making deviations more apparent and easier to interpret.

# Examples
```julia
# Compute centered L function from L function
clf = centered_l_function(l_func)

# Direct computation from data
clf = centered_l_function(data, region, radii=0.1:0.1:2.0, nk=(32,32), kmax=(0.5,0.5), tapers=tapers)
```
"""
struct CenteredLFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function CenteredLFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        return new{E, D, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end
get_short_base_estimate_name(::Type{<:CenteredLFunction}) = "centered L"
get_base_estimate_name(::Type{<:CenteredLFunction}) = "centered L function (L(r)-r)"

## required interface

computed_from(::Type{<:CenteredLFunction{E, D}}) where {E, D} = LFunction{E, D}

function allocate_estimate_memory(
        ::Type{<:CenteredLFunction}, ::Type{<:LFunction}, relevant_memory; kwargs...)
    return relevant_memory, nothing
end

extract_relevant_memory(::Type{<:CenteredLFunction}, est::LFunction) = get_estimates(est)
function extract_relevant_memory(
        ::Type{<:CenteredLFunction}, mem::EstimateMemory{<:LFunction})
    return mem.output_memory
end

validate_core_parameters(::Type{<:CenteredLFunction}; kwargs...) = nothing

function validate_memory_compatibility(
        ::Type{<:CenteredLFunction}, mem, arg::LFunction; kwargs...)
    @argcheck size(mem.output_memory) == size(get_estimates(arg))
    @argcheck eltype(mem.output_memory) == eltype(get_estimates(arg))
    return nothing
end

function resolve_missing_parameters(::Type{<:CenteredLFunction}, arg; kwargs...)
    return kwargs # no parameters beyond those in C function
end

function compute_estimate!(
        ::Type{<:CenteredLFunction{E}}, mem, source::LFunction{E}; kwargs...) where {E}
    processinfo = get_process_information(source)
    estimationinfo = get_estimation_information(source)
    radii = get_evaluation_points(source)
    value = mem.output_memory
    _l_to_centered_l_transform!(value, radii, get_estimates(source), process_trait(source))
    return CenteredLFunction{E}(radii, value, processinfo, estimationinfo)
end

get_evaluation_points(f::CenteredLFunction) = f.radii

get_estimates(f::CenteredLFunction) = f.value

## internals

"""
    _l_to_centered_l_transform!(value, radii, l_values, ::MultipleVectorTrait)

Transform L function values to centered L function for multiple vector processes.
"""
function _l_to_centered_l_transform!(
        value, radii, l_values::AbstractArray, ::MultipleVectorTrait)
    for idx in CartesianIndices(size(l_values)[1:(ndims(l_values) - 1)])
        for (i, radius) in enumerate(radii)
            value[idx, i] = _apply_centering(radius, l_values[idx, i])
        end
    end
    return value
end

"""
    _l_to_centered_l_transform!(value, radii, l_values, ::Union{MultipleTupleTrait, SingleProcessTrait})

Transform L function values to centered L function for single process or tuple traits.
"""
function _l_to_centered_l_transform!(value, radii, l_values::AbstractArray,
        ::Union{MultipleTupleTrait, SingleProcessTrait})
    for (i, radius) in enumerate(radii)
        value[i] = _apply_centering(radius, l_values[i])
    end
    return value
end

"""
    _apply_centering(radius, l_value)

Apply centering transformation: ``CL(r) = L(r) - r``.
"""
function _apply_centering(radius, l_value)
    return l_value .- radius
end
