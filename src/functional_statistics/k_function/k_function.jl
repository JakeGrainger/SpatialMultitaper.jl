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

get_evaluation_points(f::KFunction) = f.radii
get_estimates(f::KFunction) = f.value

"""
    k_function(data, region; kwargs...) -> KFunction
    k_function(data::SpatialData; kwargs...) -> KFunction
    k_function(c::CFunction; kwargs...) -> KFunction
    k_function(spectrum::NormalOrRotationalSpectra; kwargs...) -> KFunction

Compute Ripley's K function from spatial data, C functions, or spectral estimates.

Ripley's K function K(r) measures the expected number of points within distance r of a
typical point, normalized by the intensity. It is fundamental in spatial statistics for
detecting clustering or regularity in point patterns. The K function is derived from the
C function using the relationship K(r) = C(r)/λ² + V_d r^d, where λ is the process
intensity and V_d is the volume of a unit ball in d dimensions.

# Arguments
- `data`: Spatial data for direct K function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `c::CFunction`: C function estimate to transform to K function
- `spectrum::NormalOrRotationalSpectra`: Spectral estimate (converted via C function)

# Keywords
When computing directly from data, all keywords from [`c_function`](@ref) are supported:
- `radii::AbstractVector`: Distances at which to evaluate the K function (required)
- `wavenumber_radii`: Radial wavenumbers for rotational averaging
- `rotational_method`: Smoothing kernel for rotational averaging
- Plus all keywords from [`spectra`](@ref): `nk`, `kmax`, `dk`, `tapers`, `nw`, `mean_method`

# Returns
- `KFunction`: A K function estimate containing:
  - `radii`: Distance values where K function is evaluated
  - `value`: K function values (always real and positive)
  - `processinformation`: Information about the analyzed processes
  - `estimationinformation`: Details about the estimation procedure

# Mathematical Details
For a stationary point process with intensity λ, Ripley's K function is defined as:

K_{ij}(r) = λ⁻¹ E[number of additional `i` points within distance `r` of a typical `j` point]

The transformation from C function to K function is:
K(r) = C(r)/λ² + V_d r^d

where:
- C(r) is the corresponding C function value
- λ is the process intensity (mean product of the processes)
- V_d is the volume of a unit ball in d dimensions:
  - V₁ = 2 (length of unit interval)
  - V₂ = π (area of unit circle)
  - V₃ = 4π/3 (volume of unit sphere)
  - V_d = π^(d/2) / Γ(d/2 + 1) (general formula)

# Notes
- K function values are always real and non-negative
- For Poisson processes, K(r) = V_d r^d (no clustering or regularity)
- Values above V_d r^d indicate clustering; below indicate regularity
- The function automatically handles intensity normalization from process information
- Use [`partial_k_function`](@ref) for partial K functions from partial estimates

# Examples
```julia
# Direct computation from data and region
radii = 0.1:0.1:2.0
kf = k_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing C function
cf = c_function(data; radii = radii, kmax = 0.5)
kf = k_function(cf)

# From spectral estimate (via C function)
spec = spectra(data; kmax = 0.5, nw = 3)
kf = k_function(spec; radii = radii)

# Compare to theoretical Poisson process
theoretical_poisson = π .* radii.^2  # For 2D case
clustering_measure = kf.value .- theoretical_poisson
```

See also: [`c_function`](@ref), [`partial_k_function`](@ref), [`spectra`](@ref)
"""
k_function

function k_function(data, region; kwargs...)::KFunction
    return k_function(spatial_data(data, region); kwargs...)
end

function k_function(data::SpatialData; kwargs...)::KFunction
    return k_function!(c_function(data; kwargs...))
end

function k_function(est::AbstractEstimate; kwargs...)::KFunction
    mem = deepcopy(est)
    return k_function!(mem; kwargs...)
end
function k_function!(c::CFunction{E, D})::KFunction{E, D} where {E, D}
    mean_prod = get_process_information(c).mean_product
    radii = get_evaluation_points(c)
    value = _c_to_k_transform!(
        radii, get_estimates(c), process_trait(c), mean_prod, Val{D}())
    processinfo = get_process_information(c)
    estimationinfo = get_estimation_information(c)
    return KFunction{E}(radii, value, processinfo, estimationinfo)
end

function k_function!(spectrum::NormalOrRotationalSpectra; kwargs...)::KFunction
    return k_function!(c_function(spectrum; kwargs...))
end

# Partial K functions

"""
    partial_k_function(data, region; kwargs...) -> KFunction{PartialTrait}
    partial_k_function(data::SpatialData; kwargs...) -> KFunction{PartialTrait}
    partial_k_function(spectrum::NormalOrRotationalSpectra; kwargs...) -> KFunction{PartialTrait}
    partial_k_function(c::CFunction{PartialTrait}) -> KFunction{PartialTrait}

Compute partial Ripley's K function from spatial data or estimates.

The partial K function represents the direct spatial clustering or regularity patterns
after removing the linear influence of all other processes. It is computed from partial
C functions using the same transformation as regular K functions, but applied to partial
estimates that show only direct relationships between processes.

# Arguments
- `data`: Spatial data for direct partial K function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `spectrum`: Spectral estimate (partial or marginal trait)
- `c::CFunction{PartialTrait}`: Partial C function estimate

# Keywords
All keywords from [`partial_c_function`](@ref) and [`k_function`](@ref) are supported:
- `radii::AbstractVector`: Distances at which to evaluate the partial K function (required)
- Plus all spectral estimation keywords: `nk`, `kmax`, `dk`, `tapers`, `nw`, `mean_method`

# Returns
- `KFunction{PartialTrait}`: A partial K function estimate with the same structure as
    regular K functions but representing direct spatial relationships.

# Notes
- Automatically converts marginal spectra/C functions to partial when needed
- Uses the same intensity normalization and volume correction as regular K functions
- Results represent direct spatial clustering after removing indirect effects
- Cannot compute from marginal C functions (must use partial C functions or convert first)

# Examples
```julia
# Direct computation from data
radii = 0.1:0.1:2.0
partial_kf = partial_k_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing partial C function
partial_cf = partial_c_function(data; radii = radii, kmax = 0.5)
partial_kf = partial_k_function(partial_cf)

# From marginal spectrum (automatically converts to partial)
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
partial_kf = partial_k_function(marginal_spec; radii = radii)

# Compare direct vs total clustering
regular_kf = k_function(data; radii = radii, kmax = 0.5)
direct_clustering = partial_kf.value .- regular_kf.value
```

See also: [`k_function`](@ref), [`partial_c_function`](@ref), [`partial_spectra`](@ref)
"""
partial_k_function

function partial_k_function(data, region; kwargs...)::KFunction{PartialTrait}
    return partial_k_function(spatial_data(data, region); kwargs...)
end

function partial_k_function(data::SpatialData; kwargs...)::KFunction{PartialTrait}
    return k_function!(partial_c_function(data; kwargs...))
end

function partial_k_function(est::AbstractEstimate; kwargs...)::KFunction{PartialTrait}
    mem = deepcopy(est)
    return partial_k_function!(mem; kwargs...)
end

function partial_k_function!(spectrum::NormalOrRotationalSpectra{PartialTrait};
        kwargs...)::KFunction{PartialTrait}
    return k_function!(spectrum; kwargs...)
end

function partial_k_function!(spectrum::NormalOrRotationalSpectra{MarginalTrait};
        kwargs...)::KFunction{PartialTrait}
    return k_function!(partial_spectra!(spectrum); kwargs...)
end

partial_k_function!(c::CFunction{PartialTrait}) = k_function!(c)

function partial_k_function!(::CFunction{MarginalTrait})
    throw(ArgumentError(
        "Cannot compute partial K function from marginal C function. " *
        "Compute from partial spectral estimates or use partial_c_function first."
    ))
end

# Internal transformation functions

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
