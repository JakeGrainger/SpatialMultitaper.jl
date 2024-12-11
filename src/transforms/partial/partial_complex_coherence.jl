struct PartialComplexCoherence{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	partial_coherence::N
	function PartialComplexCoherence(freq, partial_coherence)
		P = checkfreqdomaininputs(freq, partial_coherence)
		new{P, typeof(freq), typeof(partial_coherence)}(freq, partial_coherence)
	end
end
getfreq(est::PartialComplexCoherence) = est.freq
getestimate(est::PartialComplexCoherence) = est.partial_coherence

function partial_complex_coherence(x::AbstractMatrix)
	-complex_coherence(inv(x))
end

function partial_complex_coherence(spectrum::SpectralEstimate)
	return PartialComplexCoherence(
		spectrum.freq,
		apply_transform(partial_complex_coherence, spectrum.power),
	)
end
