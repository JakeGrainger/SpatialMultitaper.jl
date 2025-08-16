struct PartialPhase{D,F,P,N} <: AnisotropicEstimate{D,P}
    freq::NTuple{D,F}
    partial_phase::N
    function PartialPhase(freq::NTuple{D,F}, partial_phase) where {D,F}
        P = checkinputs(freq, partial_phase)
        new{D,F,P,typeof(partial_phase)}(freq, partial_phase)
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
    )
end
