import SpatialMultitaper: centroid, grid2side, CartesianGrid, Point, Box, unitless_spacing, unitless_origin, unitless_minimum, unitless_measure
@testset "utils" begin
    @testset "unitless" begin
        grid = CartesianGrid((1.2,4.5), (1.3,20.4), dims=(10,10))
        @test unitless_spacing(grid) == ((1.3-1.2)/10, (20.4-4.5)/10)
        @test unitless_origin(grid) == (1.2, 4.5)
        @test unitless_minimum(grid) == (1.2, 4.5)
        @test unitless_measure(Box(Point(0,0), Point(1,1))) == 1
    end
    @testset "grid2side" begin
        grid = CartesianGrid((1.2,4.5), (1.3,20.4), dims=(10,10))
        gridsides = grid2side(grid)
        @test centroid(grid,1) == Point(gridsides[1][1], gridsides[2][1])
        @test centroid(grid,2) == Point(gridsides[1][2], gridsides[2][1])
        @test centroid(grid,11) == Point(gridsides[1][1], gridsides[2][2])
        @test centroid(grid,100) == Point(gridsides[1][end], gridsides[2][end])
    end
end