struct PartialMagnitudeCoherence2{D,F,P,N} <: FrequencyDomainEstimate{D,P}
    freq::NTuple{D,F}
    partial_coherence::N
    function PartialMagnitudeCoherence2(freq::NTuple{D,F}, partial_coherence) where {D,F}
        P = checkfreqdomaininputs(freq, partial_coherence)
        new{D,F,P,typeof(partial_coherence)}(freq, partial_coherence)
    end
end
getfreq(est::PartialMagnitudeCoherence2) = est.freq
getestimate(est::PartialMagnitudeCoherence2) = est.partial_coherence

function partial_magnitude_coherence2(x::AbstractMatrix)
    abs2.(partial_complex_coherence(x))
end

function partial_magnitude_coherence2(spectrum::SpectralEstimate)
    return PartialMagnitudeCoherence2(
        spectrum.freq,
        apply_transform(partial_magnitude_coherence2, spectrum.power),
    )
end
