struct PartialMagnitudeCoherence{D,F,P,N} <: FrequencyDomainEstimate{D,P}
    freq::NTuple{D,F}
    partial_coherence::N
    function PartialMagnitudeCoherence(freq::NTuple{D,F}, partial_coherence) where {D,F}
        P = checkfreqdomaininputs(freq, partial_coherence)
        new{D,F,P,typeof(partial_coherence)}(freq, partial_coherence)
    end
end
getfreq(est::PartialMagnitudeCoherence) = est.freq
getestimate(est::PartialMagnitudeCoherence) = est.partial_coherence

function partial_magnitude_coherence(x::AbstractMatrix)
    abs.(partial_complex_coherence(x))
end

function partial_magnitude_coherence(spectrum::SpectralEstimate)
    return PartialMagnitudeCoherence(
        spectrum.freq,
        apply_transform(partial_magnitude_coherence, spectrum.power),
    )
end
