struct Phase{D,F,P,N} <: FrequencyDomainEstimate{D,P}
    freq::NTuple{D,F}
    phase::N
    function Phase(freq::NTuple{D,F}, phase) where {D,F}
        P = checkfreqdomaininputs(freq, phase)
        new{D,F,P,typeof(phase)}(freq, phase)
    end
end
getfreq(est::Phase) = est.freq
getestimate(est::Phase) = est.phase

function phase(x::AbstractMatrix)
    angle.(complex_coherence(x))
end

function phase(spectrum::SpectralEstimate)
    return Phase(spectrum.freq, apply_transform(phase, spectrum.power))
end
