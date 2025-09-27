struct Coherence{E, D, P, Q, N, A, T, IP, IE} <: AbstractEstimate{E, D, P, Q, N}
    freq::A
    coherence::T
    processinformation::IP
    estimationinformation::IE
    function Coherence{E}(
            freq::NTuple{N}, coherence::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E <: EstimateTrait, D, N, T, IE}
        P, Q = checkinputs(freq, coherence, processinfo)
        IP = typeof(processinfo)
        A = typeof(freq)
        new{E, D, P, Q, N, A, T, IP, IE}(freq, coherence, processinfo, estimationinfo)
    end
    function Coherence{E}( # for inputs that are rotational
            freq::A, coherence::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E <: EstimateTrait, D, A, T, IE}
        P, Q = checkinputs(freq, coherence, processinfo)
        IP = typeof(processinfo)
        new{E, D, P, Q, 1, A, T, IP, IE}(freq, coherence, processinfo, estimationinfo)
    end
end
const RotationalCoherence{E, D, P, Q, S <: Coherence} = RotationalEstimate{E, D, P, Q, S}
getargument(est::Coherence) = est.freq
getestimate(est::Coherence) = est.coherence

function coherence(x::AbstractMatrix)
    d = diagm(sqrt.(inv.(diag(x))))
    return d * x * d
end
coherence(x::Number) = one(typeof(x))
function coherence(spectrum::Union{Spectra{E}, RotationalSpectra{E}}) where {E}
    return Coherence{E}(
        getargument(spectrum),
        apply_transform(coherence, getargument(spectrum), getestimate(spectrum)),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end

coherence(data, region; kwargs...) = coherence(spatial_data(data, region); kwargs...)
coherence(data::SpatialData; kwargs...) = coherence(spectra(data; kwargs...))

function partial_coherence(x::AbstractMatrix)
    -coherence(inv(x))
end
partial_coherence(x::Number) = one(typeof(x))
function partial_coherence(spectrum::Union{
        Spectra{E}, RotationalSpectra{E}}) where {E <: PartialTrait}
    return coherence(spectrum) # partial coherence is just coherence of the partial spectra
end

function partial_coherence(spectrum::Union{ # partial coherence from marginal spectra
        Spectra{E}, RotationalSpectra{E}}) where {E <: MarginalTrait}
    return Coherence{MarginalTrait}(
        getargument(spectrum),
        apply_transform(partial_coherence, getargument(spectrum), getestimate(spectrum)),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end
function partial_coherence(data, region; kwargs...)
    partial_coherence(spatial_data(data, region); kwargs...)
end
function partial_coherence(data::SpatialData; kwargs...)
    partial_coherence(spectra(data; kwargs...))
end

magnitude_coherence(spectrum::Spectra) = abs(coherence(spectrum))
magnitude_coherence(coh::Coherence) = abs(coh)
function magnitude_coherence(args...; kwargs...)
    magnitude_coherence(spectra(args...; kwargs...))
end

function partial_magnitude_coherence(spectrum::Spectra{MarginalTrait})
    magnitude_coherence(spectrum)
end
function partial_magnitude_coherence(spectrum::Spectra{PartialTrait})
    magnitude_coherence(partial_spectra(spectrum))
end
function partial_magnitude_coherence(coh::Coherence{MarginalTrait})
    throw(partial_from_marginal_error(Coherence, typeof(coh)))
end
partial_magnitude_coherence(coh::Coherence{PartialTrait}) = magnitude_coherence(coh)
function partial_magnitude_coherence(args...; kwargs...)
    partial_magnitude_coherence(spectra(args...; kwargs...))
end

magnitude_squared_coherence(spectrum::Spectra) = abs2(coherence(spectrum))
magnitude_squared_coherence(coh::Coherence) = abs2(coh)
function magnitude_squared_coherence(args...; kwargs...)
    magnitude_squared_coherence(spectra(args...; kwargs...))
end

function partial_magnitude_squared_coherence(spectrum::Spectra{MarginalTrait})
    magnitude_squared_coherence(spectrum)
end
function partial_magnitude_squared_coherence(spectrum::Spectra{PartialTrait})
    magnitude_squared_coherence(partial_spectra(spectrum))
end
function partial_magnitude_squared_coherence(coh::Coherence{MarginalTrait})
    throw(partial_from_marginal_error(Coherence, typeof(coh)))
end
function partial_magnitude_squared_coherence(coh::Coherence{PartialTrait})
    magnitude_squared_coherence(coh)
end
function partial_magnitude_squared_coherence(args...; kwargs...)
    partial_magnitude_squared_coherence(spectra(args...; kwargs...))
end

phase(spectrum::Union{Spectra, Coherence}) = angle(spectrum)
phase(args...; kwargs...) = phase(spectra(args...; kwargs...))

function partial_phase(spectrum::Union{Spectra{PartialTrait}, Coherence{PartialTrait}})
    phase(spectrum)
end
partial_phase(spectrum::Spectra{MarginalTrait}) = phase(partial_spectra(spectrum))
function partial_phase(coh::Coherence{MarginalTrait})
    throw(partial_from_marginal_error(Coherence, typeof(coh)))
end
partial_phase(args...; kwargs...) = partial_phase(spectra(args...; kwargs...))
