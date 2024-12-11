struct MagnitudeCoherence2{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	coherence::N
	function MagnitudeCoherence2(freq, coherence)
		P = checkfreqdomaininputs(freq, coherence)
		new{P, typeof(freq), typeof(coherence)}(freq, coherence)
	end
end
getfreq(est::MagnitudeCoherence2) = est.freq
getestimate(est::MagnitudeCoherence2) = est.coherence

function magnitude_coherence2(x::AbstractMatrix)
	abs2.(complex_coherence(x))
end

function magnitude_coherence2(spectrum::SpectralEstimate)
	return MagnitudeCoherence2(
		spectrum.freq,
		apply_transform(magnitude_coherence2, spectrum.power),
	)
end
