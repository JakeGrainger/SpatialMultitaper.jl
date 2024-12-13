struct PartialComplexCoherence{D,F,P,N} <: FrequencyDomainEstimate{D,P}
    freq::NTuple{D,F}
    partial_coherence::N
    function PartialComplexCoherence(freq::NTuple{D,F}, partial_coherence) where {D,F}
        P = checkfreqdomaininputs(freq, partial_coherence)
        new{D,F,P,typeof(partial_coherence)}(freq, partial_coherence)
    end
end
getfreq(est::PartialComplexCoherence) = est.freq
getestimate(est::PartialComplexCoherence) = est.partial_coherence

function partial_complex_coherence(x::AbstractMatrix)
    -complex_coherence(inv(x))
end

function partial_complex_coherence(spectrum::SpectralEstimate)
    return PartialComplexCoherence(
        spectrum.freq,
        apply_transform(partial_complex_coherence, spectrum.power),
    )
end
