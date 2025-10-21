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
        ::Type{KFunction}, ::Type{<:CFunction}, previous_memory; kwargs...)
    return previous_memory
end

function extract_allocation_memory end

function validate_core_parameters end

function validate_memory_compatibility end

function resolve_missing_parameters end

function apply_parameter_defaults end

function compute_estimate! end

get_evaluation_points(f::KFunction) = f.radii

get_estimates(f::KFunction) = f.value
