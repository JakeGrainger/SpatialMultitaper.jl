"""
    KFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}

Ripley's K function estimate for spatial point pattern analysis.

Ripley's K function ``K(r)`` measures the expected number of points within distance r
of a typical point, normalized by the intensity. It is fundamental in spatial
statistics for detecting clustering or regularity in point patterns.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `P`, `Q`: Process dimensions
- `A`: Type of radii array
- `T`: Type of K function values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Fields
- `radii::A`: Distances at which the K function is evaluated
- `value::T`: K function values
- `processinformation::IP`: Information about the processes
- `estimationinformation::IE`: Information about the estimation procedure

# Mathematical Background
For a stationary point process with intensity λ, Ripley's K function is:
``Kᵢⱼ(r) = λ⁻¹ E``[number of additional `i` points within distance `r` of a typical `j` point]

The relationship to the C function is:
``K(r) = C(r)/λ² + V_d r^d`` where V_d is the volume of a unit ball in d dimensions.

# Examples
```julia
# Compute K function from spectral estimate
kf = k_function(spectrum, radii=0.1:0.1:2.0)

# Direct computation from data
kf = k_function(data, region, radii=0.1:0.1:2.0, nfreq=(32,32), fmax=(0.5,0.5), tapers=tapers)
```
"""
struct KFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function KFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        P, Q = checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        return new{E, D, P, Q, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end

getbaseestimatename(::Type{<:KFunction}) = "K function"

"""
    getargument(f::KFunction)

Get the radii at which Ripley's K function is evaluated.
"""
getargument(f::KFunction) = f.radii

"""
    getestimate(f::KFunction)

Get the K function values.
"""
getestimate(f::KFunction) = f.value

"""
    k_function(data, region; kwargs...)

Compute Ripley's K function directly from data and region.
"""
function k_function(data, region; kwargs...)::KFunction
    return k_function(spatial_data(data, region); kwargs...)
end

"""
    k_function(data::SpatialData; kwargs...)

Compute Ripley's K function from spatial data.

First computes the C function, then transforms it to the K function using the
relationship ``K(r) = C(r)/λ² + V_d r^d``.

# Arguments
- `data::SpatialData`: Input spatial data
- `kwargs...`: Additional arguments passed to C function computation

# Returns
A `KFunction` object containing Ripley's K function estimates.
"""
function k_function(data::SpatialData; kwargs...)::KFunction
    return k_function(c_function(data; kwargs...))
end

"""
    k_function(c::CFunction{E, D}) where {E, D}

Compute Ripley's K function from a C function estimate.

Transforms the C function C(r) to Ripley's K function using:
``K(r) = C(r)/λ² + V_d r^d``
where ``λ`` is the process intensity and ``V_d`` is the volume of a unit d-ball.

# Arguments
- `c::CFunction`: Input C function estimate

# Returns
A `KFunction` object with the corresponding K function values.
"""
function k_function(c::CFunction{E, D})::KFunction{E, D} where {E, D}
    mean_prod = getprocessinformation(c).mean_product
    radii = getargument(c)
    value = _c_to_k_transform(radii, getestimate(c), process_trait(c), mean_prod, Val{D}())
    processinfo = getprocessinformation(c)
    estimationinfo = getestimationinformation(c)
    return KFunction{E}(radii, value, processinfo, estimationinfo)
end

"""
    k_function(spectrum::NormalOrRotationalSpectra; kwargs...)

Compute Ripley's K function from a spectral estimate.
"""
function k_function(spectrum::NormalOrRotationalSpectra; kwargs...)::KFunction
    return k_function(c_function(spectrum; kwargs...))
end

# Partial K functions

"""
    partial_k_function(data, region; kwargs...)

Compute partial Ripley's K function directly from data and region.
"""
function partial_k_function(data, region; kwargs...)::KFunction{PartialTrait}
    return partial_k_function(spatial_data(data, region); kwargs...)
end

"""
    partial_k_function(data::SpatialData; kwargs...)

Compute partial Ripley's K function from spatial data.

Partial K functions remove the linear influence of other processes.
"""
function partial_k_function(data::SpatialData; kwargs...)::KFunction{PartialTrait}
    return k_function(partial_c_function(data; kwargs...))
end

function partial_k_function(spectrum::NormalOrRotationalSpectra{PartialTrait};
        kwargs...)::KFunction{PartialTrait}
    return k_function(spectrum; kwargs...)
end

function partial_k_function(spectrum::NormalOrRotationalSpectra{MarginalTrait};
        kwargs...)::KFunction{PartialTrait}
    return k_function(partial_spectra(spectrum); kwargs...)
end

partial_k_function(c::CFunction{PartialTrait}) = k_function(c)

function partial_k_function(::CFunction{MarginalTrait})
    throw(ArgumentError(
        "Cannot compute partial K function from marginal C function. " *
        "Compute from partial spectral estimates or use partial_c_function first."
    ))
end

# Internal transformation functions

"""
    _c_to_k_transform(radii, c_values, trait, mean_prod, ::Val{D})

Transform C function values to K function values for multiple vector traits.

Handles the case where we have multiple processes stored as arrays.
"""
function _c_to_k_transform(radii, c_values::AbstractArray,
        ::MultipleVectorTrait, mean_prod, ::Val{D}) where {D}
    output = similar(c_values)

    for idx in CartesianIndices(size(c_values)[1:(ndims(c_values) - 1)])
        mean_prod_slice = mean_prod[idx]
        for (i, radius) in enumerate(radii)
            output[idx, i] = _compute_k_from_c(
                radius, c_values[idx, i], mean_prod_slice, Val{D}())
        end
    end
    return output
end

"""
    _c_to_k_transform(radii, c_values, trait, mean_prod, ::Val{D})

Transform C function values to K function values for single process or tuple traits.
"""
function _c_to_k_transform(radii, c_values::AbstractArray,
        ::Union{MultipleTupleTrait, SingleProcessTrait}, mean_prod, ::Val{D}) where {D}
    output = similar(c_values)

    for (i, radius) in enumerate(radii)
        output[i] = _compute_k_from_c(radius, c_values[i], mean_prod, Val{D}())
    end
    return output
end

"""
    _compute_k_from_c(radius, c_value, mean_prod, ::Val{D})

Core computation for transforming C to K function values.

Implements ``K(r) = C(r)/λ² + V_d r^d`` where ``V_d`` is the volume of a unit d-ball.
"""
function _compute_k_from_c(radius, c_value, mean_prod, ::Val{D}) where {D}
    unit_ball_volume = unitless_measure(Ball(Point(ntuple(_ -> 0, Val{D}())), 1))
    return c_value ./ mean_prod .+ (unit_ball_volume * (radius^D))
end

"""
    _compute_k_from_c(radius, c_value, mean_prod, ::Val{2})

Optimized computation for 2D case where unit ball volume is π.
"""
function _compute_k_from_c(radius, c_value, mean_prod, ::Val{2})
    return c_value ./ mean_prod .+ (π * radius^2)
end
