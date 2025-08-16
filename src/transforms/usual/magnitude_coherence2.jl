struct MagnitudeCoherence2{D,F,P,N} <: AnisotropicEstimate{D,P}
    freq::NTuple{D,F}
    coherence::N
    function MagnitudeCoherence2(freq::NTuple{D,F}, coherence) where {D,F}
        P = checkinputs(freq, coherence)
        new{D,F,P,typeof(coherence)}(freq, coherence)
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
    )
end
