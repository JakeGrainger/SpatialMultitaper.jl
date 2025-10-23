"""
    KFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}

Ripley's K function estimate for spatial point pattern analysis.

Ripley's K function K(r) measures the expected number of points within distance r of a
typical point, normalized by the intensity. It is fundamental in spatial statistics for
detecting clustering (K > theoretical) or regularity (K < theoretical) in point patterns.
The K function is derived from spectral estimates via the C function transformation.

# Type Parameters
- `E <: EstimateTrait`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension of the underlying process
- `A`: Type of the radii array
- `T`: Type of the K function values (typically `Vector{Float64}` or similar)
- `IP`: Type of process information structure
- `IE`: Type of estimation information structure

# Fields
- `radii::A`: Distance values at which the K function is evaluated
- `value::T`: K function values corresponding to each radius (always real and non-negative)
- `processinformation::IP`: Information about the analyzed processes (includes intensity)
- `estimationinformation::IE`: Details about the estimation procedure

# Mathematical Background
For a stationary point process with intensity λ, Ripley's K function is defined as:

K_{ij}(r) = λ⁻¹ E[number of additional `i` points within distance `r` of a typical `j` point]

The transformation from C function to K function is:
K(r) = C(r)/λ² + V_d r^d

where:
- C(r) is the corresponding C function value
- λ is the process intensity (derived from mean product)
- V_d is the volume of a unit d-dimensional ball
- For Poisson processes: K(r) = V_d r^d (baseline for comparison)

# Interpretation
- K(r) > V_d r^d: clustering at distance r
- K(r) < V_d r^d: regularity/inhibition at distance r
- K(r) = V_d r^d: random (Poisson) pattern

# Notes
- K function values are always real and non-negative
- Supports both marginal and partial variants via the trait system
- Intensity normalization is handled automatically from process information
- Common baseline comparisons: K_Poisson(r) = πr² (2D), K_Poisson(r) = (4π/3)r³ (3D)

See also: [`k_function`](@ref), [`partial_k_function`](@ref), [`CFunction`](@ref)
"""
struct KFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function KFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        return new{E, D, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end
get_short_base_estimate_name(::Type{<:KFunction}) = "K"
get_base_estimate_name(::Type{<:KFunction}) = "K function"

## required interface

computed_from(::Type{<:KFunction{E}}) where {E} = CFunction{E}

function allocate_estimate_memory(
        ::Type{<:KFunction}, ::Type{<:CFunction}, extracted_memory; kwargs...)
    return extracted_memory
end

extract_relevant_memory(::Type{KFunction}, est::CFunction) = get_estimates(est)
function extract_relevant_memory(::Type{KFunction}, mem::EstimateMemory{<:CFunction})
    return mem.estimate_type.spatial_output
end

validate_core_parameters(::Type{<:KFunction}, kwargs...) = nothing

function validate_memory_compatibility(::Type{<:KFunction}, previous_memory)
    # TODO: fill in once C function memory is done
    return nothing
end

function resolve_missing_parameters(::Type{<:KFunction}, arg; kwargs...)
    return kwargs # no parameters beyond those in C function
end

function compute_estimate!(
        ::Type{<:KFunction{E, D}}, mem, source::CFunction{E, D}; kwargs...) where {E, D}
    processinfo = get_process_information(source)
    estimationinfo = get_estimation_information(source)
    radii = get_evaluation_points(source)
    cfunc = get_estimates(source)
    mean_prod = processinfo.mean_product
    value = _c_to_k_transform!(radii, cfunc, process_trait(source), mean_prod, Val{D}())
    return KFunction{E}(radii, value, processinfo, estimationinfo)
end

get_evaluation_points(f::KFunction) = f.radii

get_estimates(f::KFunction) = f.value

## internals

"""
    _c_to_k_transform!(radii, c_values, trait, mean_prod, ::Val{D})

Transform C function values to K function values for multiple vector process traits.

Applies the transformation K(r) = C(r)/λ² + V_d r^d element-wise for multiple processes
stored as arrays, where each process pair has its own intensity λ and the volume
correction V_d r^d depends on the spatial dimension D.
"""
function _c_to_k_transform!(radii, c_values::AbstractArray,
        ::MultipleVectorTrait, mean_prod, ::Val{D}) where {D}
    for idx in CartesianIndices(size(c_values)[1:(ndims(c_values) - 1)])
        mean_prod_slice = mean_prod[idx]
        for (i, radius) in enumerate(radii)
            c_values[idx, i] = _compute_k_from_c(
                radius, c_values[idx, i], mean_prod_slice, Val{D}())
        end
    end
    return c_values
end

"""
    _c_to_k_transform!(radii, c_values, trait, mean_prod, ::Val{D})

Transform C function values to K function values for single process or tuple traits.

Applies the transformation K(r) = C(r)/λ² + V_d r^d for single processes or tuple-based
multiple processes, where the intensity λ and volume correction are applied uniformly.
"""
function _c_to_k_transform!(radii, c_values::AbstractArray,
        ::Union{MultipleTupleTrait, SingleProcessTrait}, mean_prod, ::Val{D}) where {D}
    for (i, radius) in enumerate(radii)
        c_values[i] = _compute_k_from_c(radius, c_values[i], mean_prod, Val{D}())
    end
    return c_values
end

"""
    _compute_k_from_c(radius, c_value, mean_prod, ::Val{D})

Core computation for transforming C function to K function values.

Implements the fundamental relationship K(r) = C(r)/λ² + V_d r^d where:
- C(r) is the input C function value
- λ² is the mean product (intensity squared)
- V_d is the volume of a unit d-dimensional ball
- r^d provides the appropriate scaling with radius
"""
function _compute_k_from_c(radius, c_value, mean_prod, ::Val{D}) where {D}
    unit_ball_volume = unitless_measure(Ball(Point(ntuple(_ -> 0, Val{D}())), 1))
    return c_value ./ mean_prod .+ (unit_ball_volume * (radius^D))
end

"""
    _compute_k_from_c(radius, c_value, mean_prod, ::Val{2})

Optimized computation for 2D case where the unit circle has area π.

This specialized method avoids the general unit ball volume calculation for the common
2D case, directly using π as the unit circle area in the transformation
K(r) = C(r)/λ² + π r².
"""
function _compute_k_from_c(radius, c_value, mean_prod, ::Val{2})
    return c_value ./ mean_prod .+ (π * radius^2)
end
