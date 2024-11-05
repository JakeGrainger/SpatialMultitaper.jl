function partial_covariance_density(mt_est::SpectralEstimate; filter = nothing)
    p_spectra = partial_spectra(mt_est)
    partial_cov = spectra_to_cov(p_spectra.freq, p_spectra.partial_spectra, filter)
    partial_cov_jackknifed = if isnothing(p_spectra.partial_spectra_jackknifed)
        nothing
    else
        spectra_to_cov.(Ref(p_spectra.freq), p_spectra.partial_spectra_jackknifed, Ref(filter))
    end
    lag = fftshift.(fftfreq.(length.(p_spectra.freq), inv.(step.(p_spectra.freq))))
    return (lag = lag, partial_cov = partial_cov, partial_cov_jackknifed = partial_cov_jackknifed)
end

spectra_to_cov(freq, power, ::Nothing) = (length(power) * prod(step.(freq))) .* fftshift(ifft(ifftshift(power, 3:ndims(power)), 3:ndims(power)), 3:ndims(power))
function spectra_to_cov(freq, power, filter)
    filter_values = filter.(Iterators.ProductIterator(freq))
    filtered_power = reshape(filter_values, (1, 1, size(filter_values)...)) .* power
    spectra_to_cov(freq, filtered_power, nothing)
end