"""
    CenteredLFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}

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
- `P`, `Q`: Process dimensions
- `A`: Type of radii array
- `T`: Type of centered L function values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Mathematical Background
For a Poisson process, ``L(r) = r`` theoretically. The centering transformation:
``LCL(r) = L(r) - r`` removes this trend, making deviations more apparent and easier to interpret.

# Examples
```julia
# Compute centered L function from L function
clf = centered_l_function(l_func)

# Direct computation from data
clf = centered_l_function(data, region, radii=0.1:0.1:2.0, nk=(32,32), kmax=(0.5,0.5), tapers=tapers)
```
"""
struct CenteredLFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function CenteredLFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        P, Q = checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        return new{E, D, P, Q, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end

getbaseestimatename(::Type{<:CenteredLFunction}) = "centered L function (L(r)-r)"

"""
    getargument(f::CenteredLFunction)

Get the radii at which the centered L function is evaluated.
"""
getargument(f::CenteredLFunction) = f.radii

"""
    getestimate(f::CenteredLFunction)

Get the centered L function values ``CL(r) = L(r) - r``.
"""
getestimate(f::CenteredLFunction) = f.value

"""
    centered_l_function(data, region; kwargs...)

Compute centered L function directly from data and region.
"""
function centered_l_function(data, region; kwargs...)::CenteredLFunction
    return centered_l_function(spatial_data(data, region); kwargs...)
end

"""
    centered_l_function(data::SpatialData; kwargs...)

Compute centered L function from spatial data via L function transformation.
"""
function centered_l_function(data::SpatialData; kwargs...)::CenteredLFunction
    return centered_l_function(l_function(data; kwargs...))
end

"""
    centered_l_function(spectrum::NormalOrRotationalSpectra; kwargs...)

Compute centered L function from spectral estimates.
"""
function centered_l_function(
        spectrum::NormalOrRotationalSpectra; kwargs...)::CenteredLFunction
    return centered_l_function(l_function(spectrum; kwargs...))
end

"""
    centered_l_function(c::CFunction)

Compute centered L function from C function via L function transformation.
"""
centered_l_function(c::CFunction) = centered_l_function(l_function(c))

"""
    centered_l_function(k::KFunction)

Compute centered L function from K function via L function transformation.
"""
centered_l_function(k::KFunction) = centered_l_function(l_function(k))

"""
    centered_l_function(l::LFunction{E, D}) where {E, D}

Transform L function to centered L function.

Applies the centering transformation ``CL(r) = L(r) - r`` to remove the
theoretical Poisson expectation.

# Arguments
- `l::LFunction`: Input L function estimate

# Returns
A `CenteredLFunction` object with values ``CL(r) = L(r) - r``.
"""
function centered_l_function(l::LFunction{E, D})::CenteredLFunction{E, D} where {E, D}
    radii = getargument(l)
    value = _l_to_centered_l_transform(radii, getestimate(l), process_trait(l))
    processinfo = getprocessinformation(l)
    estimationinfo = getestimationinformation(l)
    return CenteredLFunction{E}(radii, value, processinfo, estimationinfo)
end

function partial_centered_l_function(data, region; kwargs...)
    partial_centered_l_function(spatial_data(data, region); kwargs...)
end
function partial_centered_l_function(data::SpatialData; kwargs...)
    centered_l_function(partial_l_function(data; kwargs...))
end
function partial_centered_l_function(
        spectrum::NormalOrRotationalSpectra{PartialTrait}; kwargs...)
    k_function(spectrum; kwargs...)
end
function partial_centered_l_function(
        spectrum::NormalOrRotationalSpectra{MarginalTrait}; kwargs...)
    k_function(partial_spectra(spectrum); kwargs...)
end
partial_centered_l_function(est::CFunction{PartialTrait}) = centered_l_function(est)
partial_centered_l_function(est::KFunction{PartialTrait}) = centered_l_function(est)
partial_centered_l_function(est::LFunction{PartialTrait}) = centered_l_function(est)
function partial_centered_l_function(est::Union{
        LFunction{MarginalTrait}, CFunction{MarginalTrait}, KFunction{MarginalTrait}})
    throw(partial_from_marginal_error(LFunction, typeof(est)))
end

## internals

"""
    _l_to_centered_l_transform(radii, l_values, ::MultipleVectorTrait)

Transform L function values to centered L function for multiple vector processes.
"""
function _l_to_centered_l_transform(radii, l_values::AbstractArray, ::MultipleVectorTrait)
    output = zeros(eltype(l_values), size(l_values))
    for idx in CartesianIndices(size(l_values)[1:(ndims(l_values) - 1)])
        for (i, radius) in enumerate(radii)
            output[idx, i] = _apply_centering(radius, l_values[idx, i])
        end
    end
    return output
end

"""
    _l_to_centered_l_transform(radii, l_values, ::Union{MultipleTupleTrait, SingleProcessTrait})

Transform L function values to centered L function for single process or tuple traits.
"""
function _l_to_centered_l_transform(radii, l_values::AbstractArray,
        ::Union{MultipleTupleTrait, SingleProcessTrait})
    output = zeros(eltype(l_values), size(l_values))
    for (i, radius) in enumerate(radii)
        output[i] = _apply_centering(radius, l_values[i])
    end
    return output
end

"""
    _apply_centering(radius, l_value)

Apply centering transformation: ``CL(r) = L(r) - r``.
"""
function _apply_centering(radius, l_value)
    return l_value .- radius
end
