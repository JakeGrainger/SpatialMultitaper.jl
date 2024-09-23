@testset "Diagonalisation" begin
    R = Bool[
    0 0 0 0 0
    0 1 1 1 0
    0 1 1 1 0
    0 1 1 0 0
    0 0 0 0 0
    ]
    K = SpatialMultitaper.SlepianSolver.ifftshift(Bool[
        0 0 0 0 0
        0 1 1 1 0
        0 1 1 1 0
        0 1 1 1 0
        0 0 0 0 0
    ])
    h, λ = SpatialMultitaper.SlepianSolver.compute_eigenfunctions(R,K,2)
    @test h isa Vector{Matrix{Float64}}
    @test sum(abs2, h[1]) ≈ 1
    @test sum(abs2, h[2]) ≈ 1
    @test sum(x*y for (x,y) in zip(h[1],h[2])) ≈ 0 atol = 1e-5
end