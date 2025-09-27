"""
    Coherence{E, D, P, Q, N, A, T, IP, IE} <: AbstractEstimate{E, D, P, Q, N}

A coherence estimate structure containing frequency information and coherence values.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `P`, `Q`: Process dimensions
- `N`: Number of frequency dimensions
- `A`: Type of frequency argument
- `T`: Type of coherence values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Fields
- `freq`: Frequency grid or values
- `coherence`: Coherence estimates
- `processinformation`: Information about the processes
- `estimationinformation`: Information about the estimation procedure
"""
struct Coherence{E, D, P, Q, N, A, T, IP, IE} <: AbstractEstimate{E, D, P, Q, N}
    freq::A
    coherence::T
    processinformation::IP
    estimationinformation::IE
    function Coherence{E}(
            freq::NTuple{N}, coherence::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E <: EstimateTrait, D, N, T, IE}
        P, Q = checkinputs(freq, coherence, processinfo)
        IP = typeof(processinfo)
        A = typeof(freq)
        return new{E, D, P, Q, N, A, T, IP, IE}(
            freq, coherence, processinfo, estimationinfo)
    end
    function Coherence{E}( # for inputs that are rotational
            freq::A, coherence::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E <: EstimateTrait, D, A, T, IE}
        P, Q = checkinputs(freq, coherence, processinfo)
        IP = typeof(processinfo)
        return new{E, D, P, Q, 1, A, T, IP, IE}(
            freq, coherence, processinfo, estimationinfo)
    end
end
const RotationalCoherence{E, D, P, Q, S <: Coherence} = RotationalEstimate{E, D, P, Q, S}
getargument(est::Coherence) = est.freq
getestimate(est::Coherence) = est.coherence

"""
    coherence(x::AbstractMatrix)

Compute coherence from a spectral matrix by normalizing with diagonal elements.

The coherence γᵢⱼ between processes i and j is computed as:
γᵢⱼ = Sᵢⱼ / √(Sᵢᵢ * Sⱼⱼ)

# Arguments
- `x::AbstractMatrix`: Spectral matrix

# Returns
Coherence matrix with values bounded between 0 and 1 for the magnitude.
"""
function coherence(x::AbstractMatrix)
    d = diagm(sqrt.(inv.(diag(x))))
    return d * x * d
end

"""
    coherence(x::Number)

Coherence for a single process (always returns 1).
"""
coherence(x::Number) = one(typeof(x))

"""
    coherence(spectrum::NormalOrRotationalSpectra{E}) where {E}

Compute coherence from a spectral estimate.

# Arguments
- `spectrum`: A `Spectra` or `RotationalSpectra` estimate

# Returns
A `Coherence` object containing the coherence estimates.

# Throws
- `ArgumentError`: If the spectrum does not have equal input and output process sets

# Examples
```julia
# Compute coherence from spectral estimate
coh = coherence(spec)

# Direct computation from data
coh = coherence(data, region, nfreq=(32, 32), fmax=(0.5, 0.5), tapers=tapers)
```
"""
function coherence(spectrum::NormalOrRotationalSpectra{E})::Coherence{E} where {E}
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Coherence computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(spectrum))."
        ))
    end

    return _compute_coherence_estimate(spectrum, coherence, E)
end

"""
    coherence(data, region; kwargs...)

Compute coherence estimate directly from data and region.
"""
function coherence(data, region; kwargs...)::Coherence
    return coherence(spatial_data(data, region); kwargs...)
end

"""
    coherence(data::SpatialData; kwargs...)

Compute coherence estimate from spatial data.
"""
function coherence(data::SpatialData; kwargs...)::Coherence
    return coherence(spectra(data; kwargs...))
end

"""
    partial_coherence(x::AbstractMatrix)

Compute partial coherence from a spectral matrix.

Partial coherence removes the linear effects of all other processes.
"""
function partial_coherence(x::AbstractMatrix)
    return -coherence(inv(x))
end

"""
    partial_coherence(x::Number)

Partial coherence for a single process (always returns 1).
"""
partial_coherence(x::Number) = one(typeof(x))

"""
    partial_coherence(spectrum::NormalOrRotationalSpectra{PartialTrait})

Compute partial coherence from partial spectral estimates.

For partial spectral estimates, partial coherence is equivalent to coherence.
"""
function partial_coherence(spectrum::NormalOrRotationalSpectra{PartialTrait})::Coherence{PartialTrait}
    return coherence(spectrum) # partial coherence is just coherence of the partial spectra
end

"""
    partial_coherence(spectrum::NormalOrRotationalSpectra{MarginalTrait})

Compute partial coherence from marginal spectral estimates.

# Arguments
- `spectrum`: A marginal spectral estimate

# Returns
A `Coherence{PartialTrait}` object containing partial coherence estimates.

# Throws
- `ArgumentError`: If the spectrum does not have equal input and output process sets
"""
function partial_coherence(spectrum::NormalOrRotationalSpectra{MarginalTrait})::Coherence{PartialTrait}
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Partial coherence computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(spectrum))."
        ))
    end

    return _compute_coherence_estimate(spectrum, partial_coherence, PartialTrait)
end

"""
    partial_coherence(data, region; kwargs...)

Compute partial coherence estimate directly from data and region.
"""
function partial_coherence(data, region; kwargs...)::Coherence{PartialTrait}
    return partial_coherence(spatial_data(data, region); kwargs...)
end

"""
    partial_coherence(data::SpatialData; kwargs...)

Compute partial coherence estimate from spatial data.
"""
function partial_coherence(data::SpatialData; kwargs...)::Coherence{PartialTrait}
    return partial_coherence(spectra(data; kwargs...))
end

"""
    magnitude_coherence(spectrum::Spectra)

Compute magnitude coherence |γᵢⱼ| from a spectral estimate.
"""
magnitude_coherence(spectrum::Spectra) = abs(coherence(spectrum))

"""
    magnitude_coherence(coh::Coherence)

Compute magnitude of existing coherence estimate.
"""
magnitude_coherence(coh::Coherence) = abs(coh)

"""
    magnitude_coherence(args...; kwargs...)

Compute magnitude coherence directly from data arguments.
"""
function magnitude_coherence(args...; kwargs...)
    return magnitude_coherence(spectra(args...; kwargs...))
end

"""
    partial_magnitude_coherence(spectrum::Spectra{MarginalTrait})

Compute partial magnitude coherence from marginal spectral estimates.
"""
function partial_magnitude_coherence(spectrum::Spectra{MarginalTrait})
    return magnitude_coherence(spectrum)
end

function partial_magnitude_coherence(spectrum::Spectra{PartialTrait})
    return magnitude_coherence(partial_spectra(spectrum))
end

function partial_magnitude_coherence(coh::Coherence{MarginalTrait})
    throw(_partial_from_marginal_error("partial magnitude coherence", typeof(coh)))
end

partial_magnitude_coherence(coh::Coherence{PartialTrait}) = magnitude_coherence(coh)

function partial_magnitude_coherence(args...; kwargs...)
    return partial_magnitude_coherence(spectra(args...; kwargs...))
end

"""
    magnitude_squared_coherence(spectrum::Spectra)

Compute magnitude-squared coherence |γᵢⱼ|² from a spectral estimate.

Also known as the coherency or coherence function.
"""
magnitude_squared_coherence(spectrum::Spectra) = abs2(coherence(spectrum))

magnitude_squared_coherence(coh::Coherence) = abs2(coh)

function magnitude_squared_coherence(args...; kwargs...)
    return magnitude_squared_coherence(spectra(args...; kwargs...))
end

function partial_magnitude_squared_coherence(spectrum::Spectra{MarginalTrait})
    return magnitude_squared_coherence(spectrum)
end

function partial_magnitude_squared_coherence(spectrum::Spectra{PartialTrait})
    return magnitude_squared_coherence(partial_spectra(spectrum))
end

function partial_magnitude_squared_coherence(coh::Coherence{MarginalTrait})
    throw(_partial_from_marginal_error("partial magnitude squared coherence", typeof(coh)))
end

function partial_magnitude_squared_coherence(coh::Coherence{PartialTrait})
    return magnitude_squared_coherence(coh)
end

function partial_magnitude_squared_coherence(args...; kwargs...)
    return partial_magnitude_squared_coherence(spectra(args...; kwargs...))
end

"""
    phase(spectrum::Union{Spectra, Coherence})

Compute phase of coherence or cross-spectral estimates.

Returns the phase angle in radians.
"""
phase(spectrum::Union{Spectra, Coherence}) = angle(spectrum)

phase(args...; kwargs...) = phase(spectra(args...; kwargs...))

function partial_phase(spectrum::Union{Spectra{PartialTrait}, Coherence{PartialTrait}})
    return phase(spectrum)
end

partial_phase(spectrum::Spectra{MarginalTrait}) = phase(partial_spectra(spectrum))

function partial_phase(coh::Coherence{MarginalTrait})
    throw(_partial_from_marginal_error("partial phase", typeof(coh)))
end

partial_phase(args...; kwargs...) = partial_phase(spectra(args...; kwargs...))

# Helper functions

"""
    _compute_coherence_estimate(spectrum, transform_func, trait_type)

Internal helper to compute coherence estimates with proper error handling.
"""
function _compute_coherence_estimate(spectrum, transform_func, trait_type)
    trait = process_trait(spectrum)
    est = getestimate(spectrum)
    process_info = getprocessinformation(spectrum)
    estimation_info = getestimationinformation(spectrum)
    transformed = apply_transform(transform_func, est, trait)
    return Coherence{trait_type}(
        getargument(spectrum), transformed, process_info, estimation_info)
end

"""
    _partial_from_marginal_error(operation_name, estimate_type)

Generate a descriptive error for unsupported partial operations on marginal estimates.
"""
function _partial_from_marginal_error(operation_name::String, estimate_type::Type)
    return ArgumentError(
        "Cannot compute $(operation_name) from a marginal coherence estimate. " *
        "Compute from spectral estimates or use partial spectral estimates instead. " *
        "Got estimate of type: $(estimate_type)"
    )
end
