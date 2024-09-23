function transform_spectral_estimate(transform::T, mt_est::SpectralEstimate{F, P, Nothing}) where {T, F, P}
	transformed_power = mapslices(transform, mt_est.power, dims = (1, 2))
	return (freq=mt_est.freq, transformed_power=transformed_power, transformed_power_jackknifed=nothing)
end

function transform_spectral_estimate(transform::T, mt_est::SpectralEstimate{F, P, Vector{P}}) where {T, F, P}
	transformed_power = mapslices(transform, mt_est.power, dims = (1, 2))
	transformed_power_jackknifed = mapslices.(Ref(transform), mt_est.power_jackknifed, dims = (1, 2))
	return (freq=mt_est.freq, transformed_power=transformed_power, transformed_power_jackknifed=transformed_power_jackknifed)
end

macro spectraltransform(x)
    esc(quote
        $x
        function $(x.args[1].args[1])(mt_est::SpectralEstimate) 
			transformed_estimate = transform_spectral_estimate($(x.args[1].args[1]), mt_est)
			return (freq=transformed_estimate.freq, $(x.args[1].args[1])=transformed_estimate.transformed_power, $(Symbol("$(x.args[1].args[1])_jackknifed"))=transformed_estimate.transformed_power_jackknifed)
		end
    end)
end

##

@spectraltransform function complex_coherence(x::Matrix)
	return [x[i, j] / sqrt(x[i, i] * x[j, j]) for i ∈ axes(x, 1), j ∈ axes(x, 2)]
end

@spectraltransform function magnitude_coherence(x::Matrix)
	return abs.(complex_coherence(x))
end

@spectraltransform function magnitude_sq_coherence(x::Matrix)
	return abs2.(complex_coherence(x))
end

@spectraltransform function group_delay(x::Matrix)
	return angle.(complex_coherence(x))
end

@spectraltransform function partial_complex_coherence(x::Matrix)
	return -complex_coherence(inv(x))
end

@spectraltransform function partial_magnitude_coherence(x::Matrix)
	return abs.(partial_complex_coherence(x))
end

@spectraltransform function partial_magnitude_sq_coherence(x::Matrix)
	return abs2.(partial_complex_coherence(x))
end

@spectraltransform function partial_group_delay(x::Matrix)
	return angle.(partial_complex_coherence(x))
end

@spectraltransform function partial_spectra(x::Matrix)
	C = inv(x)
	par_coh = -complex_coherence(C)
	Sa = inv.(diag(C))
	return [
		i == j ? Sa[i] : par_coh[i, j] / (1 - abs2(par_coh[i, j])) * sqrt(Sa[i] * Sa[j])
		for i in axes(C, 1), j in axes(C, 2)]
end

