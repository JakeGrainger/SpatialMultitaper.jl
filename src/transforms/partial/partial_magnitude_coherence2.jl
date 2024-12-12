struct PartialMagnitudeCoherence2{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	partial_coherence::N
	function PartialMagnitudeCoherence2(freq, partial_coherence)
		P = checkfreqdomaininputs(freq, partial_coherence)
		new{P, typeof(freq), typeof(partial_coherence)}(freq, partial_coherence)
	end
end
getfreq(est::PartialMagnitudeCoherence2) = est.freq
getestimate(est::PartialMagnitudeCoherence2) = est.partial_coherence

function partial_magnitude_coherence2(x::AbstractMatrix)
	abs2.(partial_complex_coherence(x))
end

function partial_magnitude_coherence2(spectrum::SpectralEstimate)
	return PartialMagnitudeCoherence2(
		spectrum.freq,
		apply_transform(partial_magnitude_coherence2, spectrum.power),
	)
end
