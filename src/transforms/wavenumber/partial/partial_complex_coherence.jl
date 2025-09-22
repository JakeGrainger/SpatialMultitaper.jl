struct PartialComplexCoherence{F, N, I, T, D, P, Q} <: AnisotropicEstimate{D, P, Q}
    freq::NTuple{D, F}
    partial_coherence::N
    processinformation::I
    estimationinformation::T
    function PartialComplexCoherence(
            freq::NTuple{D, F}, partial_coherence, processinfo, estimationinfo) where {D, F}
        P, Q = checkinputs(freq, partial_coherence, processinfo)
        new{F, typeof(partial_coherence), typeof(processinfo),
            typeof(estimationinfo), D, P, Q}(
            freq, partial_coherence, processinfo, estimationinfo)
    end
end
getargument(est::PartialComplexCoherence) = est.freq
getestimate(est::PartialComplexCoherence) = est.partial_coherence

function partial_complex_coherence(x::AbstractMatrix)
    -complex_coherence(inv(x))
end

function partial_complex_coherence(spectrum::SpectralEstimate)
    return PartialComplexCoherence(
        spectrum.freq,
        apply_transform(partial_complex_coherence, spectrum.power),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end
