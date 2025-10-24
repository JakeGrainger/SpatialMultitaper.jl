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
