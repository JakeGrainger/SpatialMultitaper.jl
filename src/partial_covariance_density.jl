function partial_covariance_density(mt_est::SpectralEstimate)
    p_spectra = partial_spectra(mt_est)
    partial_cov = spectra_to_cov(p_spectra.freq, p_spectra.partial_spectra)
    partial_cov_jackknifed = if isnothing(p_spectra.partial_spectra_jackknifed)
        nothing
    else
        spectra_to_cov.(Ref(p_spectra.freq), p_spectra.partial_spectra_jackknifed)
    end
    lag = fftshift.(fftfreq.(length.(p_spectra.freq), inv.(step.(p_spectra.freq))))
    return (lag = lag, partial_cov = partial_cov, partial_cov_jackknifed = partial_cov_jackknifed)
end

spectra_to_cov(freq, power) = (length(power) * prod(step.(freq))) .* fftshift(ifft(ifftshift(power, 3:ndims(power)), 3:ndims(power)), 3:ndims(power))