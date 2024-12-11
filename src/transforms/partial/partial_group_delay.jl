struct PartialGroupDelay{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	partial_group_delay::N
	function PartialGroupDelay(freq, partial_group_delay)
		P = checkfreqdomaininputs(freq, partial_group_delay)
		new{P, typeof(freq), typeof(partial_group_delay)}(freq, partial_group_delay)
	end
end
getfreq(est::PartialGroupDelay) = est.freq
getestimate(est::PartialGroupDelay) = est.partial_group_delay

function partial_group_delay(x::AbstractMatrix)
	angle.(partial_complex_coherence(x))
end

function partial_group_delay(spectrum::SpectralEstimate)
	return GroupDelay(spectrum.freq, apply_transform(partial_group_delay, spectrum.power))
end
