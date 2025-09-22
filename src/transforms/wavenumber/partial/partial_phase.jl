struct PartialPhase{F, N, I, T, D, P, Q} <: AnisotropicEstimate{D, P, Q}
    freq::NTuple{D, F}
    partial_phase::N
    processinformation::I
    estimationinformation::T
    function PartialPhase(
            freq::NTuple{D, F}, partial_phase, processinfo, estimationinfo) where {D, F}
        P, Q = checkinputs(freq, partial_phase, processinfo)
        new{F, typeof(partial_phase), typeof(processinfo),
            typeof(estimationinfo), D, P, Q}(
            freq, partial_phase, processinfo, estimationinfo)
    end
end
getargument(est::PartialPhase) = est.freq
getestimate(est::PartialPhase) = est.partial_phase

function partial_phase(x::AbstractMatrix)
    angle.(partial_complex_coherence(x))
end

function partial_phase(spectrum::SpectralEstimate)
    return PartialPhase(
        spectrum.freq,
        apply_transform(partial_phase, spectrum.power),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end
