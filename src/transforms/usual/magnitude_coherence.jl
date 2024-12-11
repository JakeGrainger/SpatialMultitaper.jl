struct MagnitudeCoherence{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	coherence::N
	function MagnitudeCoherence(freq, coherence)
		P = checkfreqdomaininputs(freq, coherence)
		new{P, typeof(freq), typeof(coherence)}(freq, coherence)
	end
end
getfreq(est::MagnitudeCoherence) = est.freq
getestimate(est::MagnitudeCoherence) = est.coherence

function magnitude_coherence(x::AbstractMatrix)
	abs.(complex_coherence(x))
end

function magnitude_coherence(spectrum::SpectralEstimate)
	return MagnitudeCoherence(
		spectrum.freq,
		apply_transform(magnitude_coherence, spectrum.power),
	)
end
