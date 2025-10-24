"""
    LFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}

L function estimate derived from Ripley's K function for spatial analysis.

The L function is a variance-stabilizing transformation of Ripley's K function defined as:
L(r) = sign(K(r)) × (|K(r)|/V_d)^(1/d), where V_d is the volume of a unit ball in d
dimensions. This transformation linearizes the expected behavior for Poisson processes
(L(r) = r) and provides better statistical properties for hypothesis testing and pattern
detection compared to the K function.

# Type Parameters
- `E <: EstimateTrait`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension of the underlying process
- `A`: Type of the radii array
- `T`: Type of the L function values (typically `Vector{Float64}` or similar)
- `IP`: Type of process information structure
- `IE`: Type of estimation information structure

# Fields
- `radii::A`: Distance values at which the L function is evaluated
- `value::T`: L function values corresponding to each radius (real-valued, can be negative)
- `processinformation::IP`: Information about the analyzed processes
- `estimationinformation::IE`: Details about the estimation procedure

# Mathematical Background
The L function transformation provides several advantages over K functions:

**Variance Stabilization**: The power transformation (1/d) stabilizes the variance across
different scales, making statistical tests more reliable.

**Linear Baseline**: For Poisson processes, L(r) = r provides a simple linear baseline
for comparison, unlike K(r) = V_d r^d which is hard to compare visually.

# Statistical Interpretation
- **L(r) = r**: Random (Poisson-like) pattern at distance r
- **L(r) > r**: Clustering at distance r (more points than expected)
- **L(r) < r**: Regularity/inhibition at distance r (fewer points than expected)

# Notes
- Better variance properties than K functions across scales
- Supports both marginal and partial variants via the trait system
- Commonly used baseline: L(r) = r for any dimension (vs K baseline V_d r^d)
- Particularly useful for statistical testing and confidence envelope construction

See also: [`l_function`](@ref), [`partial_l_function`](@ref), [`LFunction`](@ref)
"""
struct LFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function LFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        return new{E, D, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end

get_short_base_estimate_name(::Type{<:LFunction}) = "L"
get_base_estimate_name(::Type{<:LFunction}) = "L function"

## required interface

computed_from(::Type{<:LFunction{E}}) where {E} = KFunction{E}

function allocate_estimate_memory(
        ::Type{<:LFunction}, ::Type{<:KFunction}, relevant_memory; kwargs...)
    return relevant_memory
end

extract_relevant_memory(::Type{LFunction}, est::KFunction) = get_estimates(est)
function extract_relevant_memory(::Type{LFunction}, mem::EstimateMemory{<:KFunction})
    return mem.output_memory
end

validate_core_parameters(::Type{<:LFunction}, kwargs...) = nothing

function validate_memory_compatibility(::Type{<:LFunction}, mem, arg::KFunction; kwargs...)
    @argcheck size(mem.output_memory) == size(get_estimates(arg))
    @argcheck eltype(mem.output_memory) == eltype(get_estimates(arg))
    return nothing
end

function resolve_missing_parameters(::Type{<:LFunction}, arg; kwargs...)
    return kwargs # no parameters beyond those in C function
end

function compute_estimate!(
        ::Type{<:LFunction{E, D}}, mem, source::KFunction{E, D}; kwargs...) where {E, D}
    processinfo = get_process_information(source)
    estimationinfo = get_estimation_information(source)
    radii = get_evaluation_points(source)
    value = _k_to_l_transform!(mem.output_memory, get_estimates(source), Val{D}())
    return LFunction{E}(radii, value, processinfo, estimationinfo)
end

get_evaluation_points(f::LFunction) = f.radii

get_estimates(f::LFunction) = f.value

## internals

"""
    _k_to_l_transform!(value, k_values::AbstractArray, ::Val{D})

Transform K function values to L function values for any array structure.

Applies the transformation L(r) = sign(K) × (|K|/V_d)^(1/d) element-wise to the entire
array, preserving the array structure while converting each K value to its corresponding
L value with appropriate sign handling for negative K values.
"""
function _k_to_l_transform!(value, k_values::AbstractArray, ::Val{D}) where {D}
    for idx in eachindex(k_values)
        value[idx] = _compute_l_from_k(k_values[idx], Val{D}())
    end
    return value
end

"""
    _compute_l_from_k(k_value, ::Val{D})

Core L function computation with sign preservation: L(r) = sign(K) × (|K|/V_d)^(1/d).

This transformation provides variance stabilization while preserving the sign of K values
to handle cases where strong regularity produces negative K values. The absolute value
ensures the power operation is well-defined, while the sign multiplication preserves
the directional information about clustering vs regularity.

# Arguments
- `k_value`: K function value(s) to transform
- `::Val{D}`: Spatial dimension for unit ball volume calculation
"""
function _compute_l_from_k(k_value, ::Val{D}) where {D}
    unit_ball_volume = unitless_measure(Ball(Point(ntuple(_ -> 0, Val{D}())), 1))
    return sign.(k_value) .* (abs.(k_value) ./ unit_ball_volume) .^ (1 / D)
end

"""
    _compute_l_from_k(k_value, ::Val{2})

Optimized L function computation for 2D: L(r) = sign(K) × √(|K|/π).

This specialized method avoids the general unit ball volume calculation for the common
2D case, directly using π as the unit circle area and the square root operation
(power 1/2) for the variance-stabilizing transformation in 2D spatial analysis.
"""
function _compute_l_from_k(k_value, ::Val{2})
    return sign.(k_value) .* sqrt.(abs.(k_value) ./ π)
end
