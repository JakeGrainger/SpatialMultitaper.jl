struct ComplexCoherence{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	coherence::N
	function ComplexCoherence(freq, coherence)
		P = checkfreqdomaininputs(freq, coherence)
		new{P, typeof(freq), typeof(coherence)}(freq, coherence)
	end
end
getfreq(est::ComplexCoherence) = est.freq
getestimate(est::ComplexCoherence) = est.coherence

function complex_coherence(x::AbstractMatrix)
	d = diagm(sqrt.(inv.(diag(x))))
	return d * x * d
end

function complex_coherence(spectrum::SpectralEstimate)
	return ComplexCoherence(
		spectrum.freq,
		apply_transform(complex_coherence, spectrum.power),
	)
end
