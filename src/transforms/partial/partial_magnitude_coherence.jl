struct PartialMagnitudeCoherence{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	partial_coherence::N
	function PartialMagnitudeCoherence(freq, partial_coherence)
		P = checkfreqdomaininputs(freq, partial_coherence)
		new{P, typeof(freq), typeof(partial_coherence)}(freq, partial_coherence)
	end
end
getfreq(est::PartialMagnitudeCoherence) = est.freq
getestimate(est::PartialMagnitudeCoherence) = est.coherence

function partial_magnitude_coherence(x::AbstractMatrix)
	abs.(partial_complex_coherence(x))
end

function partial_magnitude_coherence(spectrum::SpectralEstimate)
	return PartialMagnitudeCoherence(
		spectrum.freq,
		apply_transform(partial_magnitude_coherence, spectrum.power),
	)
end
