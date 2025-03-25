struct MagnitudeCoherence{D,F,P,N} <: FrequencyDomainEstimate{D,P}
    freq::NTuple{D,F}
    coherence::N
    function MagnitudeCoherence(freq::NTuple{D,F}, coherence) where {D,F}
        P = checkfreqdomaininputs(freq, coherence)
        new{D,F,P,typeof(coherence)}(freq, coherence)
    end
end
getfreq(est::MagnitudeCoherence) = est.freq
getestimate(est::MagnitudeCoherence) = est.coherence

function magnitude_coherence(x::AbstractMatrix)
    abs.(complex_coherence(x))
end

function magnitude_coherence(spectrum::SpectralEstimate)
    return MagnitudeCoherence(
        spectrum.freq,
        apply_transform(magnitude_coherence, spectrum.power),
    )
end
