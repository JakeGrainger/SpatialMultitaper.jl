struct GroupDelay{D,F,P,N} <: FrequencyDomainEstimate{D,P}
    freq::NTuple{D,F}
    group_delay::N
    function GroupDelay(freq::NTuple{D,F}, group_delay) where {D,F}
        P = checkfreqdomaininputs(freq, group_delay)
        new{D,F,P,typeof(group_delay)}(freq, group_delay)
    end
end
getfreq(est::GroupDelay) = est.freq
getestimate(est::GroupDelay) = est.group_delay

function group_delay(x::AbstractMatrix)
    angle.(complex_coherence(x))
end

function group_delay(spectrum::SpectralEstimate)
    return GroupDelay(spectrum.freq, apply_transform(group_delay, spectrum.power))
end
