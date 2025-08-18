pspec = partial_spectra(mt_est)
pspec_alt = partial_spectra(mt_est, 1, 1, 2:3, 2:3)
@test pspec[1, 1].partial_spectra ≈ pspec_alt.partial_spectra

function long_partial_spectra(x)
    C = inv(x)
    
    return [
        i == j ? 1 / C[i,i] : -C[i,j] / (C[i,i] * C[j,j]) for
        i in axes(C, 1), j in axes(C, 2)
    ]
end
@test long_partial_spectra(mt_est.power[1, 1]) ≈ pspec.partial_spectra[1, 1]