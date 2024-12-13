struct PartialGroupDelay{D,F,P,N} <: FrequencyDomainEstimate{D,P}
    freq::NTuple{D,F}
    partial_group_delay::N
    function PartialGroupDelay(freq::NTuple{D,F}, partial_group_delay) where {D,F}
        P = checkfreqdomaininputs(freq, partial_group_delay)
        new{D,F,P,typeof(partial_group_delay)}(freq, partial_group_delay)
    end
end
getfreq(est::PartialGroupDelay) = est.freq
getestimate(est::PartialGroupDelay) = est.partial_group_delay

function partial_group_delay(x::AbstractMatrix)
    angle.(partial_complex_coherence(x))
end

function partial_group_delay(spectrum::SpectralEstimate)
    return PartialGroupDelay(
        spectrum.freq,
        apply_transform(partial_group_delay, spectrum.power),
    )
end
