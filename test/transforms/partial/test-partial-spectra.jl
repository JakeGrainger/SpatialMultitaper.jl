mt_est = make_simple_example()
pspec = Spmt.partial_spectra_uncorrected(mt_est)

function long_partial_spectra(x)
    C = inv(x)
    
    return [
        i == j ? 1 / C[i,i] : -C[i,j] / (C[i,i] * C[j,j] - abs2(C[i,j])) for
        i in axes(C, 1), j in axes(C, 2)
    ]
end

function full_partial_spectra(x)
    return [
        x[i,j] - (x[i:i,Spmt.Not(i,j)] / x[Spmt.Not(i,j), Spmt.Not(i,j)] * x[Spmt.Not(i,j), j:j])[1] for
        i in axes(x, 1), j in axes(x, 2)
    ]
end

@test long_partial_spectra(mt_est.power[1, 1]) ≈ pspec.partial_spectra[1, 1] # 1,1 means first wavenumber
@test full_partial_spectra(mt_est.power[1, 1]) ≈ pspec.partial_spectra[1, 1]
@test partial_spectra(mt_est.power[1, 1][1:2,1:2], nothing) ≈ partial_spectra(collect(mt_est.power[1, 1][1:2,1:2]), nothing) # check is specialised methods are equivalent