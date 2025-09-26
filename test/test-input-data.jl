using SpatialMultitaper, Test
import SpatialMultitaper: check_observations

@testset "check_observations" begin
    p1 = georef((rf = ones(2),), PointSet([Point(0), Point(1)]))
    p2 = georef((rf = ones(2),), PointSet([Point(0, 0), Point(1, 1)]))
    p3 = PointSet([Point(0, 0), Point(1, 1)])
    region = Box((0, 0), (1, 1))

    @test_throws ArgumentError check_observations((p1, p2), region)
    @test_throws ArgumentError check_observations((p1, p3), region)

    pattern = georef((marks = [1, 2],), PointSet([Point(0), Point(1)]))
    griddata = georef((rf = rand(6 * 6),), CartesianGrid((0, 0), (3, 3), dims = (6, 6)))
    multi_grid = georef(
        (rf1 = griddata.rf, rf2 = rand(6 * 6)),
        CartesianGrid((0, 0), (3, 3), dims = (6, 6))
    )

    @test_throws ArgumentError check_observations((pattern, griddata), region)
    @test_throws ArgumentError check_observations((domain(pattern), griddata), region)
    @test_throws ArgumentError check_observations([domain(pattern), griddata], region)
end
