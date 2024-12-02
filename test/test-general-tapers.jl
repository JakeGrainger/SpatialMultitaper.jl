@testset "downsample_spacing" begin
	g = CartesianGrid((100, 100), (1.0, 2.0), (0.1, 0.3))
	@test SpatialMultitaper.downsample_spacing(g, 4) == (0.1, 0.3) .* 4
	@test SpatialMultitaper.downsample_spacing(g, nothing) == (0.1, 0.3)
end
@testset "reprocess" begin
	g = CartesianGrid((50, 50), (1.0, 2.0), (0.1, 0.3))
	x = rand(50, 50)
	@test SpatialMultitaper.reprocess(x, g, nothing) == x
	@test downsample(SpatialMultitaper.reprocess(x, g, 10), 10) ≈ x[1:5, 1:5] ./ 10
end
@testset "pixelate_region" begin
	g = CartesianGrid((0.0, 0.0), (1.0, 1.0), dims = (100, 100))
	b = Box(Point(0, 0), Point(1, 1))
	@test SpatialMultitaper.pixelate_region(g, b) == ones(Bool, 100, 100)
	@test SpatialMultitaper.pixelate_region(g, discretize(b)) == ones(Bool, 100, 100)
	@test SpatialMultitaper.pixelate_region(grid2side(g), b) == ones(Bool, 100, 100)
end
@testset "optimaltapers" begin
	b = Box(Point(0, 0), Point(1, 1))
	g = CartesianGrid((0.0, 0.0), (1.0, 1.0), dims = (100, 100))
	g2 = CartesianGrid((0.0, 0.0), (1.0, 0.5), dims = (100, 100))
	h = SpatialMultitaper.optimaltapers(b, g, freq_region = Ball(Point(0, 0), 0.1),
		ntapers = 2,
		freq_res = 200, freq_downsample = 2, tol = 0.1)
	@test h[1] isa Vector{Matrix{Float64}}
	@test_throws AssertionError SpatialMultitaper.optimaltapers(b, g2,
		freq_region = Ball(Point(0, 0), 0.1), ntapers = 2,
		freq_res = 200, freq_downsample = 2, tol = 0.1)
end
