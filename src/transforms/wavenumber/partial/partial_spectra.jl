struct PartialSpectra{F, N, I, T, D, P, Q} <: AnisotropicEstimate{D, P, Q}
    freq::NTuple{D, F}
    partial_spectra::N
    processinformation::I
    estimationinformation::T
    function PartialSpectra(
            freq::NTuple{D, F}, partial_coherence, processinfo, estimationinfo) where {D, F}
        P, Q = checkinputs(freq, partial_coherence, processinfo)

        new{F, typeof(partial_coherence), typeof(processinfo),
            typeof(estimationinfo), D, P, Q}(
            freq, partial_coherence, processinfo, estimationinfo)
    end
end
getargument(est::PartialSpectra) = est.freq
getestimate(est::PartialSpectra) = est.partial_spectra

"""
	partial_spectra(x::Matrix)

If only a matrix `x` is passed, computes the partial spectra for all indices.
The diagonal elements are the spectra of residuals given all other processes
The i j th element (for i ≠ j) is the partial spectra between the ith and jth process given all other processes not including i and j.

	partial_spectra(x::Matrix, i1::Int, i2::Int, c1, c2)

If specific indices are requested, computes the partial spectra for the i1th index conditioned on the indices in c1 vs the i2th index conditioned on the indices in c2.
"""
function partial_spectra(x::SMatrix, ::Nothing)
    g = inv(x)
    A = diagm((diag(g)))
    g2 = abs2.(g)
    denom = A * ones(typeof(x)) * A - g2 + diagm(diag(g2))
    return (g ./ denom) .* (2I - ones(typeof(x)))
    # computes -gⱼₖ / (gⱼⱼ gₖₖ - |gⱼₖ|²) if j ≠ k
    # computes  1 / gⱼⱼ if j = k
end

function partial_spectra(x::AbstractMatrix, ::Nothing)
    C = inv(x)

    return [i == j ? 1 / C[i, i] : -C[i, j] / (C[i, i] * C[j, j] - abs2(C[i, j]))
            for
            i in axes(C, 1), j in axes(C, 2)]
end

function partial_spectra(x::SMatrix{2, 2, T, 4}, ::Nothing) where {T}
    g = inv(x)
    p11 = 1 / g[1, 1]
    p22 = 1 / g[2, 2]
    p12 = -g[1, 2] / (g[1, 1] * g[2, 2] - abs2(g[1, 2]))
    p21 = conj(p12)

    return SMatrix{2, 2, T, 4}(
        p11,
        p21,
        p12,
        p22 # column major
    )
end

function partial_spectra(x::SMatrix{Q, Q, T, N}, ntapers::Int) where {Q, T, N}
    p = partial_spectra(x, nothing)
    denom = ntapers .* ones(typeof(x)) .- Q .+ 2 - I # so that M - Q + 2 off diag and M - Q + 1 on diag
    return ntapers ./ denom .* p
end

function partial_spectra(x::AbstractMatrix{T}, ntapers::Int) where {T}
    Q = size(x, 1)
    p = partial_spectra(x, nothing)
    for i in axes(p, 1), j in axes(p, 2)
        @inbounds p[i, j] *= ntapers / (ntapers - Q + 2 - (i == j))
    end
    return p
end

function partial_spectra(spectrum::SpectralEstimate)
    return PartialSpectra(
        spectrum.freq,
        apply_transform(
            partial_spectra, spectrum.power, getestimationinformation(spectrum).ntapers),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end

function partial_spectra_uncorrected(spectrum::SpectralEstimate)
    return PartialSpectra(
        spectrum.freq,
        apply_transform(partial_spectra, spectrum.power, nothing),
        getprocessinformation(spectrum), getestimationinformation(spectrum)
    )
end
