struct Phase{F, N, I, T, D, P, Q} <: AnisotropicEstimate{D, P, Q}
    freq::NTuple{D, F}
    phase::N
    processinformation::I
    estimationinformation::T
    function Phase(freq::NTuple{D, F}, phase, processinfo, estimationinfo) where {D, F}
        P, Q = checkinputs(freq, phase, processinfo)
        new{F, typeof(phase), typeof(processinfo), typeof(estimationinfo), D, P, Q}(
            freq, phase, processinfo, estimationinfo)
    end
end
getargument(est::Phase) = est.freq
getestimate(est::Phase) = est.phase

function phase(x::AbstractMatrix)
    angle.(complex_coherence(x))
end

function phase(spectrum::SpectralEstimate)
    return Phase(spectrum.freq, apply_transform(phase, spectrum.power),
        getprocessinformation(spectrum), getestimationinformation(spectrum))
end
