@testset "downsample_spacing" begin
    g = CartesianGrid((100, 100), (1.0, 2.0), (0.1, 0.3))
    @test Spmt.downsample_spacing(g, 4) == (0.1, 0.3) .* 4
    @test Spmt.downsample_spacing(g, nothing) == (0.1, 0.3)
end

@testset "reprocess" begin
    g = CartesianGrid((50, 50), (1.0, 2.0), (0.1, 0.3))
    x = rand(50, 50)
    @test Spmt.reprocess(x, g, nothing) == x
    @test downsample(Spmt.reprocess(x, g, 10), 10) ≈ x[1:5, 1:5] ./ 10
end

@testset "pixelate_region" begin
    g = CartesianGrid((0.0, 0.0), (1.0, 1.0), dims = (100, 100))
    b = Box(Point(0, 0), Point(1, 1))

    result = zeros(Bool, 100, 100) # because the edge points will be set to zero by not including them
    result[2:(end - 1), 2:(end - 1)] .= true

    @test Spmt.pixelate_region(g, b) == result
end

@testset "optimaltapers" begin
    b = Box(Point(0, 0), Point(1, 1))
    g = CartesianGrid((0.0, 0.0), (1.0, 1.0), dims = (100, 100))
    g2 = CartesianGrid((0.0, 0.0), (1.0, 0.5), dims = (100, 100))
    options = (
        freq_region = Ball(Point(0, 0), 0.1),
        ntapers = 2,
        freq_res = 200,
        freq_downsample = 2,
        tol = 0.1
    )

    h = Spmt.optimaltapers(b, g; options...)
    @test h[1] isa Vector{Matrix{Float64}}
    @test_throws AssertionError Spmt.optimaltapers(b, g2; options...)
end
