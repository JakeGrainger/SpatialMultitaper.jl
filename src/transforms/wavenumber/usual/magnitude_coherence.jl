struct MagnitudeCoherence{F, N, I, T, D, P, Q} <: AnisotropicEstimate{D, P, Q}
    freq::NTuple{D, F}
    coherence::N
    processinformation::I
    estimationinformation::T
    function MagnitudeCoherence(
            freq::NTuple{D, F}, coherence, processinfo, estimationinfo) where {D, F}
        P, Q = checkinputs(freq, coherence, processinfo)
        new{F, typeof(coherence), typeof(processinfo), typeof(estimationinfo), D, P, Q}(
            freq, coherence, processinfo, estimationinfo)
    end
end
getargument(est::MagnitudeCoherence) = est.freq
getestimate(est::MagnitudeCoherence) = est.coherence

function magnitude_coherence(x::AbstractMatrix)
    abs.(complex_coherence(x))
end

function magnitude_coherence(spectrum::SpectralEstimate)
    return MagnitudeCoherence(
        spectrum.freq,
        apply_transform(magnitude_coherence, spectrum.power),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end
