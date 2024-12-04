@testset "unitless" begin
	grid = CartesianGrid((1.2, 4.5), (1.3, 20.4), dims = (10, 10))
	@test all(SpatialMultitaper.unitless_spacing(grid) .≈ (0.01, 1.59)) # (1.3 - 1.2)/10, (20.4 - 4.5)/10
	@test SpatialMultitaper.unitless_origin(grid) == (1.2, 4.5)
	@test SpatialMultitaper.unitless_minimum(grid) == (1.2, 4.5)
	@test SpatialMultitaper.unitless_measure(Box(Point(0, 0), Point(1, 1))) == 1
end
@testset "grid2side" begin
	grid = CartesianGrid((1.2, 4.5), (1.3, 20.4), dims = (10, 10))
	gridsides = SpatialMultitaper.grid2side(grid)
	@test centroid(grid, 1) == Point(gridsides[1][1], gridsides[2][1])
	@test centroid(grid, 2) == Point(gridsides[1][2], gridsides[2][1])
	@test centroid(grid, 11) == Point(gridsides[1][1], gridsides[2][2])
	@test centroid(grid, 100) == Point(gridsides[1][end], gridsides[2][end])
end

@testset "downsample" begin
	x = 1:10
	@test downsample(x, 2) == 1:2:9
	@test downsample(x, 2) !== 2:2:10
	@test_throws AssertionError downsample(x, 11)
	@test_throws AssertionError downsample(x, 0)
	y = rand(10, 10, 10)
	@test downsample(y, 2) == y[1:2:end, 1:2:end, 1:2:end]
	@test downsample(y, nothing) === y
end
@testset "pad" begin
	M = ones(10, 10)
	@test pad(M, 100) == [[ones(10, 10); zeros(90, 10)] zeros(100, 90)]
	@test pad(1:4, 100) == [1:4; zeros(96)]
	@test_throws AssertionError pad(M, (5, 11))
	@test SpatialMultitaper.pad(zeros(4, 4, 4), (10, 10, 10)) == zeros(10, 10, 10)
	y = SpatialMultitaper.pad(ones(2, 2, 2), (10, 10, 10))
	@test all(x -> x == 1, y[1:2, 1:2, 1:2])
	@test all(x -> x == 0, y[3:10, 3:10, 3:10])
	@test_throws AssertionError SpatialMultitaper.pad(ones(3, 5), (2, 2))
	@test SpatialMultitaper.pad(1:3, (5,)) == [1, 2, 3, 0, 0]
end
@testset "grid2side" begin
	g = CartesianGrid((100, 100), (1.0, 2.0), (0.1, 0.3))
	grid = grid2side(g)
	@test length(grid[1]) == 100
	@test length(grid[2]) == 100
	@test grid[1][1] == 1.05
	@test grid[2][1] == 2.15
	@test step.(grid) == (0.1, 0.3)
end
@testset "upsample" begin
	@testset "1d" begin
		g = CartesianGrid((100,), (1.0,), (0.1,))
		x = rand(10)
		@test SpatialMultitaper.upsample(x, g, nothing) == x
		@test_throws AssertionError SpatialMultitaper.upsample(x, g, 2)
		@test_throws AssertionError SpatialMultitaper.upsample(x, g, (2,))
		@test size(SpatialMultitaper.upsample(x, g, 10)) == (100,)
		@test downsample(SpatialMultitaper.upsample(x, g, 10), 10) ≈ x
	end
	@testset "2d" begin
		g = CartesianGrid((100, 100), (1.0, 2.0), (0.1, 0.3))
		x = rand(10, 10)
		@test SpatialMultitaper.upsample(x, g, nothing) == x
		@test_throws AssertionError SpatialMultitaper.upsample(x, g, 2)
		@test_throws AssertionError SpatialMultitaper.upsample(x, g, (2, 10))
		@test size(SpatialMultitaper.upsample(x, g, 10)) == (100, 100)
		@test downsample(SpatialMultitaper.upsample(x, g, 10), 10) ≈ x
	end
end
