struct PartialSpectra{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	partial_spectra::N
	function PartialSpectra(freq, partial_coherence)
		P = checkfreqdomaininputs(freq, partial_coherence)
		new{P, typeof(freq), typeof(partial_coherence)}(freq, partial_coherence)
	end
end
getfreq(est::PartialSpectra) = est.freq
getestimate(est::PartialSpectra) = est.partial_coherence


"""
	partial_spectra(x::Matrix)

If only a matrix `x` is passed, computes the partial spectra for all indices. 
The diagonal elements are the spectra of residuals given all other processes
The i j th element (for i ≠ j) is the partial coherence between the ith and jth process given all other processes not including i and j.
	
	partial_spectra(x::Matrix, i1::Int, i2::Int, c1, c2)

If specific indices are requested, computes the partial spectra for the i1th index conditioned on the indices in c1 vs the i2th index conditioned on the indices in c2.
"""
function partial_spectra(x::AbstractMatrix)
	C = inv(x)
	par_coh = -complex_coherence(C)
	invCd = inv.(diag(C))
	invCd_root = diagm(sqrt.(invCd))
	part1 = (invCd_root * par_coh * invCd_root)
	(part1 - diagm(diag(part1)) + diagm(invCd)) ./
	(1.0 .- abs2.(par_coh) + I)
end

function partial_spectra(x::AbstractMatrix, i1::Int, i2::Int, c1, c2) # note, this can be made more efficient
	x[i1, i2] - transpose(x[i1, c1]) / x[c1, c1] * x[c1, i2] -
	transpose(x[i1, c2]) / x[c2, c2] * x[c2, i2] +
	transpose(x[i1, c1]) / x[c1, c1] * x[c1, c2] / x[c2, c2] * x[c2, i2]
end


function partial_spectra(spectrum::SpectralEstimate)
	return PartialSpectra(spectrum.freq, apply_transform(partial_spectra, spectrum.power))
end

function partial_spectra(spectrum::SpectralEstimate, i1::Int, i2::Int, c1, c2)
	function wrapped_partial_spectra(x)
		partial_spectra(x, i1, i2, c1, c2)
	end # maybe better to make a version of apply_transform for passing args...
	return PartialSpectra(
		spectrum.freq,
		apply_transform(wrapped_partial_spectra, spectrum.power),
	)
end
