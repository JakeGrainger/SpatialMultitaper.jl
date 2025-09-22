struct PartialMagnitudeCoherence{F, N, I, T, D, P, Q} <: AnisotropicEstimate{D, P, Q}
    freq::NTuple{D, F}
    partial_coherence::N
    processinformation::I
    estimationinformation::T
    function PartialMagnitudeCoherence(
            freq::NTuple{D, F}, partial_coherence, processinfo, estimationinfo) where {D, F}
        P, Q = checkinputs(freq, partial_coherence, processinfo)
        new{F, typeof(partial_coherence), typeof(processinfo),
            typeof(estimationinfo), D, P, Q}(
            freq, partial_coherence, processinfo, estimationinfo)
    end
end
getargument(est::PartialMagnitudeCoherence) = est.freq
getestimate(est::PartialMagnitudeCoherence) = est.partial_coherence

function partial_magnitude_coherence(x::AbstractMatrix)
    abs.(partial_complex_coherence(x))
end

function partial_magnitude_coherence(spectrum::SpectralEstimate)
    return PartialMagnitudeCoherence(
        spectrum.freq,
        apply_transform(partial_magnitude_coherence, spectrum.power),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end
