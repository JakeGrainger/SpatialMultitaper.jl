struct PartialMagnitudeCoherence2{F, N, I, T, D, P, Q} <: AnisotropicEstimate{D, P, Q}
    freq::NTuple{D, F}
    partial_coherence::N
    processinformation::I
    estimationinformation::T
    function PartialMagnitudeCoherence2(
            freq::NTuple{D, F}, partial_coherence, processinfo, estimationinfo) where {D, F}
        P, Q = checkinputs(freq, partial_coherence, processinfo)
        new{F, typeof(partial_coherence), typeof(processinfo),
            typeof(estimationinfo), D, P, Q}(
            freq, partial_coherence, processinfo, estimationinfo)
    end
end
getargument(est::PartialMagnitudeCoherence2) = est.freq
getestimate(est::PartialMagnitudeCoherence2) = est.partial_coherence

function partial_magnitude_coherence2(x::AbstractMatrix)
    abs2.(partial_complex_coherence(x))
end

function partial_magnitude_coherence2(spectrum::SpectralEstimate)
    return PartialMagnitudeCoherence2(
        spectrum.freq,
        apply_transform(partial_magnitude_coherence2, spectrum.power),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end
