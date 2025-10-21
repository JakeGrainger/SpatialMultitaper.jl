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

See also: [`l_function`](@ref), [`partial_l_function`](@ref), [`KFunction`](@ref)
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
getshortbaseestimatename(::Type{<:LFunction}) = "L"
getbaseestimatename(::Type{<:LFunction}) = "L function"

getevaluationpoints(f::LFunction) = f.radii
getestimates(f::LFunction) = f.value

"""
    l_function(data, region; kwargs...) -> LFunction
    l_function(data::SpatialData; kwargs...) -> LFunction
    l_function(k::KFunction; kwargs...) -> LFunction
    l_function(est::AbstractEstimate; kwargs...) -> LFunction

Compute L function from spatial data, K functions, or other estimates.

The L function is a variance-stabilizing transformation of Ripley's K function defined as:
L(r) = (K(r)/V_d)^(1/d), where V_d is the volume of a unit ball in d dimensions. This
transformation makes the L function approximately linear for Poisson processes (L(r) ≈ r),
which greatly simplifies the interpretation of clustering and regularity patterns compared
to the K function.

# Arguments
- `data`: Spatial data for direct L function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `k::KFunction`: K function estimate to transform to L function
- `est::AbstractEstimate`: Any estimate that can be converted via K function

# Keywords
When computing directly from data, all keywords from [`k_function`](@ref) are supported:
- `radii::AbstractVector`: Distances at which to evaluate the L function (required)
- Plus all keywords from [`c_function`](@ref) and [`spectra`](@ref): `wavenumber_radii`,
    `rotational_method`, `nk`, `kmax`, `dk`, `tapers`, `nw`, `mean_method`

# Returns
- `LFunction`: An L function estimate containing:
  - `radii`: Distance values where L function is evaluated
  - `value`: L function values (real-valued, can be negative)
  - `processinformation`: Information about the analyzed processes
  - `estimationinformation`: Details about the estimation procedure

# Mathematical Details
The L function transformation from K function is:

L(r) = sign(K(r)) × (|K(r)|/V_d)^(1/d)

where:
- K(r) is the input K function value
- V_d is the volume of a unit d-dimensional ball:
  - V₁ = 2, V₂ = π, V₃ = 4π/3, V_d = π^(d/2)/Γ(d/2 + 1)
- The sign is preserved to handle negative K values properly
- The power 1/d provides variance stabilization

**Key Properties:**
- For Poisson processes: L(r) = r (linear baseline)
- L(r) > r indicates clustering at distance r
- L(r) < r indicates regularity/inhibition at distance r
- Better variance properties than K function for statistical testing

# Notes
- L function provides better variance stabilization than K function
- Linear behavior L(r) = r for Poisson processes simplifies interpretation
- Use [`partial_l_function`](@ref) for partial L functions from partial estimates

# Examples
```julia
# Direct computation from data and region
radii = 0.1:0.1:2.0
lf = l_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing K function
kf = k_function(data; radii = radii, kmax = 0.5)
lf = l_function(kf)

# From spectral estimate (via K function)
spec = spectra(data; kmax = 0.5, nw = 3)
lf = l_function(spec; radii = radii)
```

See also: [`k_function`](@ref), [`partial_l_function`](@ref), [`c_function`](@ref)
"""
function l_function(data, region; kwargs...)::LFunction
    return l_function(spatial_data(data, region); kwargs...)
end

function l_function(data::SpatialData; kwargs...)::LFunction
    return l_function!(k_function(data; kwargs...))
end

function l_function(est::AbstractEstimate; kwargs...)::LFunction
    mem = deepcopy(est)
    return l_function!(mem; kwargs...)
end

function l_function!(est::AbstractEstimate; kwargs...)::LFunction
    return l_function!(k_function!(est; kwargs...))
end

function l_function!(k::KFunction{E, D})::LFunction{E, D} where {E, D}
    radii = getevaluationpoints(k)
    value = _k_to_l_transform!(getestimates(k), Val{D}())
    processinfo = getprocessinformation(k)
    estimationinfo = getestimationinformation(k)
    return LFunction{E}(radii, value, processinfo, estimationinfo)
end

# Partial L functions

"""
    partial_l_function(data, region; kwargs...) -> LFunction{PartialTrait}
    partial_l_function(data::SpatialData; kwargs...) -> LFunction{PartialTrait}
    partial_l_function(est::AbstractEstimate; kwargs...) -> LFunction{PartialTrait}
    partial_l_function(k::KFunction{PartialTrait}) -> LFunction{PartialTrait}

Compute partial L function from spatial data or estimates.

The partial L function represents the direct spatial clustering or regularity patterns
after removing the linear influence of all other processes. It applies the same
variance-stabilizing transformation L(r) = (K(r)/V_d)^(1/d) to partial K functions,
providing a linearized view of direct spatial relationships with the same baseline
L(r) = r for direct Poisson-like behavior.

# Arguments
- `data`: Spatial data for direct partial L function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `est::AbstractEstimate`: Any estimate that can be converted via partial K function
- `k::KFunction{PartialTrait}`: Partial K function estimate to transform

# Keywords
All keywords from [`partial_k_function`](@ref) and [`l_function`](@ref) are supported:
- `radii::AbstractVector`: Distances at which to evaluate the partial L function (required)
- Plus all spectral estimation keywords: `nk`, `kmax`, `dk`, `tapers`, `nw`, `mean_method`

# Returns
- `LFunction{PartialTrait}`: A partial L function estimate with the same structure as
    regular L functions but representing direct spatial relationships.

# Mathematical Details
The partial L function transformation is identical to the regular L function:

L_partial(r) = sign(K_partial(r)) × (|K_partial(r)|/V_d)^(1/d)

but applied to partial K functions K_partial(r) that show only direct relationships.

**Interpretation:**
- L_partial(r) = r: no direct clustering or regularity at distance r
- L_partial(r) > r: direct clustering after removing indirect effects
- L_partial(r) < r: direct regularity/inhibition after removing indirect effects

# Notes
- Automatically converts marginal estimates to partial when needed
- Uses the same variance-stabilization properties as regular L functions
- Cannot compute from marginal K or C functions (must use partial estimates)
- Baseline comparison is still L(r) = r regardless of partial vs marginal

# Examples
```julia
# Direct computation from data
radii = 0.1:0.1:2.0
partial_lf = partial_l_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing partial K function
partial_kf = partial_k_function(data; radii = radii, kmax = 0.5)
partial_lf = partial_l_function(partial_kf)

# From marginal spectrum (automatically converts to partial)
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
partial_lf = partial_l_function(marginal_spec; radii = radii)

```

See also: [`l_function`](@ref), [`partial_k_function`](@ref), [`partial_spectra`](@ref)
"""
function partial_l_function(data, region; kwargs...)::LFunction{PartialTrait}
    return partial_l_function(spatial_data(data, region); kwargs...)
end

function partial_l_function(data::SpatialData; kwargs...)::LFunction{PartialTrait}
    return l_function!(partial_k_function(data; kwargs...))
end

function partial_l_function(est::AbstractEstimate; kwargs...)::LFunction{PartialTrait}
    mem = deepcopy(est)
    return partial_l_function!(mem; kwargs...)
end

function partial_l_function!(spectrum::NormalOrRotationalSpectra{PartialTrait};
        kwargs...)::LFunction{PartialTrait}
    return l_function!(k_function!(spectrum; kwargs...))
end

function partial_l_function!(spectrum::NormalOrRotationalSpectra{MarginalTrait};
        kwargs...)::LFunction{PartialTrait}
    return l_function!(k_function!(partial_spectra!(spectrum); kwargs...))
end

partial_l_function!(est::CFunction{PartialTrait}) = l_function!(est)
partial_l_function!(est::KFunction{PartialTrait}) = l_function!(est)

function partial_l_function!(est::Union{CFunction{MarginalTrait}, KFunction{MarginalTrait}})
    throw(ArgumentError(
        "Cannot compute partial L function from marginal $(typeof(est)). " *
        "Use partial spectral estimates or partial_k_function first."
    ))
end

# Internal transformation functions

"""
    _k_to_l_transform!(k_values::AbstractArray, ::Val{D})

Transform K function values to L function values for any array structure.

Applies the transformation L(r) = sign(K) × (|K|/V_d)^(1/d) element-wise to the entire
array, preserving the array structure while converting each K value to its corresponding
L value with appropriate sign handling for negative K values.
"""
function _k_to_l_transform!(k_values::AbstractArray, ::Val{D}) where {D}
    for idx in eachindex(k_values)
        k_values[idx] = _compute_l_from_k(k_values[idx], Val{D}())
    end
    return k_values
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
