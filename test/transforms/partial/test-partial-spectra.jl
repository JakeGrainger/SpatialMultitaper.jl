mt_est = make_simple_example()
pspec = Spmt.partial_spectra_uncorrected(mt_est)

function long_partial_spectra(x)
    C = inv(x)

    return [
        i == j ? 1 / C[i, i] : -C[i, j] / (C[i, i] * C[j, j] - abs2(C[i, j])) for
        i in axes(C, 1), j in axes(C, 2)
    ]
end

function full_partial_spectra(x)
    return [
        x[i, j] -
        (x[i:i, Spmt.Not(i, j)]/x[Spmt.Not(i, j), Spmt.Not(i, j)]*x[Spmt.Not(i, j), j:j])[1]
        for i in axes(x, 1), j in axes(x, 2)
    ]
end

S = mt_est.power[1, 1]
@test long_partial_spectra(S) ≈ pspec.partial_spectra[1, 1] # 1,1 means first wavenumber
@test full_partial_spectra(S) ≈ pspec.partial_spectra[1, 1]
@test partial_spectra(S[1:2, 1:2], nothing) ≈ partial_spectra(collect(S[1:2, 1:2]), nothing) # check is specialised methods are equivalent

A = Spmt.@SMatrix randn(ComplexF64, 10, 10)
M = A * A'

for n = 4:2:10
    submatrix_smat = M[Spmt.SOneTo(n), Spmt.SOneTo(n)]
    submatrix_mat = collect(submatrix_smat)
    @test submatrix_smat isa Spmt.SMatrix{n,n,ComplexF64,n^2}
    @test submatrix_mat isa Matrix{ComplexF64}
    @test submatrix_smat ≈ submatrix_mat

    p = n ÷ 2
    p_from_smat = Spmt.split_partial_spectra(submatrix_smat, nothing)
    p_from_mat = Spmt.split_partial_spectra(submatrix_mat, nothing)
    @test p_from_smat isa Spmt.SMatrix{p,p,ComplexF64,p^2}
    @test p_from_mat isa Matrix{ComplexF64}
    @test p_from_smat ≈ p_from_mat
end

@test partial_spectra(one(Spmt.SMatrix{3,3,ComplexF64,9}), 5) ≈
      (5 / 3) .* one(Spmt.SMatrix{3,3,ComplexF64,9})
@test partial_spectra(Spmt.diagm(ones(3)), 5) ≈ (5 / 3) .* Spmt.diagm(ones(3))