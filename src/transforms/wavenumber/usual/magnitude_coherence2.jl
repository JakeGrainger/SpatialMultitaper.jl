struct MagnitudeCoherence2{F, N, I, T, D, P, Q} <: AnisotropicEstimate{D, P, Q}
    freq::NTuple{D, F}
    coherence::N
    processinformation::I
    estimationinformation::T
    function MagnitudeCoherence2(
            freq::NTuple{D, F}, coherence, processinfo, estimationinfo) where {D, F}
        P, Q = checkinputs(freq, coherence, processinfo)
        new{F, typeof(coherence), typeof(processinfo), typeof(estimationinfo), D, P, Q}(
            freq, coherence, processinfo, estimationinfo)
    end
end
getargument(est::MagnitudeCoherence2) = est.freq
getestimate(est::MagnitudeCoherence2) = est.coherence

function magnitude_coherence2(x::AbstractMatrix)
    abs2.(complex_coherence(x))
end

function magnitude_coherence2(spectrum::SpectralEstimate)
    return MagnitudeCoherence2(
        spectrum.freq,
        apply_transform(magnitude_coherence2, spectrum.power),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end
