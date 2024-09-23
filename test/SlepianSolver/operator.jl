@testset "Operators" begin
    R = [
    0 0 0 0 0
    0 1 1 1 0
    0 1 1 1 0
    0 1 1 0 0
    0 0 0 0 0
    ]
    L = [
        0 0 0 0 0
        0 1 1 1 0
        0 1 1 1 0
        0 1 1 1 0
        0 0 0 0 0
    ]

    C,invV = SpatialMultitaper.SlepianSolver.make_concentration_operator(R,L,Val{false}())
    @test eltype(C) == ComplexF64
    @test size(C) == (length(R), length(L))
    x = ones(eltype(C), size(C,2))
    cx = C * x
    @test cx isa Vector{ComplexF64}
    @test length(cx) == (size(C,1))
    y = rand(eltype(C), size(C,2))
    a = 2.3
    @test C * (x.+a.*y) ≈ (C*x) .+ (a .* (C*y))

    @testset "1d case" begin
        C,invV = SpatialMultitaper.SlepianSolver.make_concentration_operator(ones(10),[zeros(4);ones(2);zeros(4)],Val{true}())
        @test C * ones(10) isa Vector{Float64}
    end

end

@testset "inv Vec" begin
    x = rand(20)
    y = rand(ComplexF64, 10, 2)
    z = rand(5,4)
    @test SpatialMultitaper.SlepianSolver.fast_reshape!(y,x) == reshape(x,size(y))
    @test SpatialMultitaper.SlepianSolver.fast_reshape!(z,x) == reshape(x,size(z))
    @test SpatialMultitaper.SlepianSolver.fast_reshape!(x,z) == vec(z) == reshape(z,size(x))
end
