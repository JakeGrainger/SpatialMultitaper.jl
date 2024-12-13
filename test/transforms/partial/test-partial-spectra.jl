pspec = partial_spectra(mt_est)
pspec_alt = partial_spectra(mt_est, 1, 1, 2:3, 2:3)
@test pspec[1, 1].partial_spectra ≈ pspec_alt.partial_spectra

function long_partial_spectra(x)
    C = inv(x)
    par_coh = -complex_coherence(C)
    Sa = inv.(Spmt.diag(C))
    return [
        i == j ? Sa[i] : par_coh[i, j] / (1 - abs2(par_coh[i, j])) * sqrt(Sa[i] * Sa[j]) for
        i in axes(C, 1), j in axes(C, 2)
    ]
end
@test long_partial_spectra(mt_est.power[1, 1]) ≈ pspec.partial_spectra[1, 1]