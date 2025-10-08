"""
    LFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}

L function estimate derived from Ripley's K function for spatial analysis.

The L function is a variance-stabilizing transformation of Ripley's K function:
``L(r) = (K(r)/V_d)^(1/d)``
where V_d is the volume of a unit ball in d dimensions.

For a Poisson process, ``L(r) = r``, making deviations easier to interpret.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `A`: Type of radii array
- `T`: Type of L function values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Mathematical Background
The L function transformation provides:
- Better variance stabilization than K function
- Linear behavior ``L(r) = r`` for Poisson processes
- Easier interpretation of clustering/regularity patterns

# Examples
```julia
# Compute L function from K function
lf = l_function(k_func)

# Direct computation from data
lf = l_function(data, region, radii=0.1:0.1:2.0, nk=(32,32), kmax=(0.5,0.5), tapers=tapers)
```
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

"""
    getargument(f::LFunction)

Get the radii at which the L function is evaluated.
"""
getargument(f::LFunction) = f.radii

"""
    getestimate(f::LFunction)

Get the L function values.
"""
getestimate(f::LFunction) = f.value

"""
    l_function(data, region; kwargs...)

Compute L function directly from data and region.
"""
function l_function(data, region; kwargs...)::LFunction
    return l_function(spatial_data(data, region); kwargs...)
end

"""
    l_function(data::SpatialData; kwargs...)

Compute L function from spatial data.
"""
function l_function(data::SpatialData; kwargs...)::LFunction
    return l_function!(k_function(data; kwargs...))
end

"""
    l_function(est::AbstractEstimate)

Compute L function from some estimate.
"""
function l_function(est::AbstractEstimate; kwargs...)::LFunction
    mem = deepcopy(est)
    return l_function!(mem; kwargs...)
end

"""
    l_function!(spec::NormalOrRotationalSpectra; kwargs...)

Compute L function from spectral estimates.
"""
function l_function!(est::AbstractEstimate; kwargs...)::LFunction
    return l_function!(k_function!(est; kwargs...))
end

"""
    l_function!(k::KFunction{E, D}) where {E, D}

Transform Ripley's K function to L function.
"""
function l_function!(k::KFunction{E, D})::LFunction{E, D} where {E, D}
    radii = getargument(k)
    value = _k_to_l_transform!(getestimate(k), Val{D}())
    processinfo = getprocessinformation(k)
    estimationinfo = getestimationinformation(k)
    return LFunction{E}(radii, value, processinfo, estimationinfo)
end

# Partial L functions

"""
    partial_l_function(data, region; kwargs...)

Compute partial L function directly from data and region.
"""
function partial_l_function(data, region; kwargs...)::LFunction{PartialTrait}
    return partial_l_function(spatial_data(data, region); kwargs...)
end

"""
    partial_l_function(data::SpatialData; kwargs...)

Compute partial L function from spatial data.
"""
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
"""
function _k_to_l_transform!(k_values::AbstractArray, ::Val{D}) where {D}
    for idx in eachindex(k_values)
        k_values[idx] = _compute_l_from_k(k_values[idx], Val{D}())
    end
    return k_values
end

"""
    _compute_l_from_k(k_value, ::Val{D})

Core L function computation: ``L(r) = (K(r)/V_d)^(1/d)``.

Handles sign preservation for negative K values by using:
``L(r) = sign(K) * (|K|/V_d)^(1/d)``
"""
function _compute_l_from_k(k_value, ::Val{D}) where {D}
    unit_ball_volume = unitless_measure(Ball(Point(ntuple(_ -> 0, Val{D}())), 1))
    return sign.(k_value) .* (abs.(k_value) ./ unit_ball_volume) .^ (1 / D)
end

"""
    _compute_l_from_k(k_value, ::Val{2})

Optimized L function computation for 2D: ``L(r) = sign(K) * sqrt(|K|/π)``.
"""
function _compute_l_from_k(k_value, ::Val{2})
    return sign.(k_value) .* sqrt.(abs.(k_value) ./ π)
end
