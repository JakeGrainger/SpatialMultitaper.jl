struct PartialSpectra{D,F,P,N} <: AnisotropicEstimate{D,P}
    freq::NTuple{D,F}
    partial_spectra::N
    function PartialSpectra(freq::NTuple{D,F}, partial_coherence) where {D,F}
        P = checkinputs(freq, partial_coherence)
        new{D,F,P,typeof(partial_coherence)}(freq, partial_coherence)
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
function partial_spectra(x::AbstractMatrix, ::Nothing)
    g = inv(x)
    A = diagm(inv.(diag(g)))
    return (A * g * A) .* (2I - ones(typeof(x)))
    # computes -gⱼₖ / (gⱼⱼ gₖₖ) if j ≠ k
    # computes  gⱼₖ / (gⱼⱼ gₖₖ) if j = k
end

function partial_spectra(x::AbstractMatrix, ntapers)
    p = partial_spectra(x, nothing)
    Q = size(x,1)
    denom = ntapers - Q + 2 - I # so that M - Q + 2 off diag and M - Q + 1 on diag
    return ntapers ./ denom .* p
end

function partial_spectra(x::AbstractMatrix, i1::Int, i2::Int, c1, c2, ::Nothing) # note, this can be made more efficient
    x[i1, i2] - transpose(x[i1, c1]) / x[c1, c1] * x[c1, i2] -
    transpose(x[i1, c2]) / x[c2, c2] * x[c2, i2] +
    transpose(x[i1, c1]) / x[c1, c1] * x[c1, c2] / x[c2, c2] * x[c2, i2]
end


function partial_spectra(spectrum::SpectralEstimate)
    return PartialSpectra(spectrum.freq, apply_transform(partial_spectra, spectrum.power, spectrum.ntapers))
end

function partial_spectra(spectrum::SpectralEstimate, i1::Int, i2::Int, c1, c2)
    function wrapped_partial_spectra(x)
        partial_spectra(x, i1, i2, c1, c2, spectrum.ntapers)
    end # maybe better to make a version of apply_transform for passing args...
    return PartialSpectra(
        spectrum.freq,
        apply_transform(wrapped_partial_spectra, spectrum.power),
    )
end

function partial_spectra_uncorrected(spectrum::SpectralEstimate)
    return PartialSpectra(spectrum.freq, apply_transform(partial_spectra_uncorrected, spectrum.power, nothing))
end

function partial_spectra_uncorrected(spectrum::SpectralEstimate, i1::Int, i2::Int, c1, c2)
    function wrapped_partial_spectra(x)
        partial_spectra(x, i1, i2, c1, c2, nothing)
    end # maybe better to make a version of apply_transform for passing args...
    return PartialSpectra(
        spectrum.freq,
        apply_transform(wrapped_partial_spectra, spectrum.power),
    )
end
