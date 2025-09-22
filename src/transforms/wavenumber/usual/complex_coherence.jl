struct ComplexCoherence{F, N, I, T, D, P, Q} <: AnisotropicEstimate{D, P, Q}
    freq::NTuple{D, F}
    coherence::N
    processinformation::I
    estimationinformation::T
    function ComplexCoherence(
            freq::NTuple{D, F}, coherence, processinfo, estimationinfo) where {D, F}
        P, Q = checkinputs(freq, coherence, processinfo)
        new{F, typeof(coherence), typeof(processinfo), typeof(estimationinfo), D, P, Q}(
            freq, coherence, processinfo, estimationinfo)
    end
end
getargument(est::ComplexCoherence) = est.freq
getestimate(est::ComplexCoherence) = est.coherence

function complex_coherence(x::AbstractMatrix)
    d = diagm(sqrt.(inv.(diag(x))))
    return d * x * d
end

function complex_coherence(spectrum::SpectralEstimate)
    return ComplexCoherence(
        spectrum.freq,
        apply_transform(complex_coherence, spectrum.power),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end
