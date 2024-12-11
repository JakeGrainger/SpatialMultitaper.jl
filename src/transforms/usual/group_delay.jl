struct GroupDelay{P, F, N} <: FrequencyDomainEstimate{P}
	freq::F
	group_delay::N
	function GroupDelay(freq, group_delay)
		P = checkfreqdomaininputs(freq, group_delay)
		new{P, typeof(freq), typeof(group_delay)}(freq, group_delay)
	end
end
getfreq(est::GroupDelay) = est.freq
getestimate(est::GroupDelay) = est.group_delay

function group_delay(x::AbstractMatrix)
	angle.(complex_coherence(x))
end

function group_delay(spectrum::SpectralEstimate)
	return GroupDelay(spectrum.freq, apply_transform(group_delay, spectrum.power))
end
