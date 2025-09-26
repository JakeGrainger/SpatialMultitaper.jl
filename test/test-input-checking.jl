using SpatialMultitaper, Test
import SpatialMultitaper: check_spatial_data

@testset "check_spatial_data" begin
    p1 = georef((rf = ones(2),), PointSet([Point(0), Point(1)]))
    p2 = georef((rf = ones(2),), PointSet([Point(0, 0), Point(1, 1)]))
    p3 = PointSet([Point(0, 0), Point(1, 1)])
    @test ((p2, p3), 2) == check_spatial_data((p2, p3))
    @test_throws AssertionError check_spatial_data((p1, p2))
    @test_throws AssertionError check_spatial_data((p1, p3))

    pattern = georef((marks = [1, 2],), PointSet([Point(0), Point(1)]))
    griddata = georef((rf = rand(6 * 6),), CartesianGrid((0, 0), (3, 3), dims = (6, 6)))
    multi_grid = georef(
        (rf1 = griddata.rf, rf2 = rand(6 * 6)),
        CartesianGrid((0, 0), (3, 3), dims = (6, 6))
    )

    @test_throws AssertionError check_spatial_data((pattern, griddata))
    @test_throws AssertionError check_spatial_data((domain(pattern), griddata))
    @test_throws AssertionError check_spatial_data([domain(pattern), griddata])
    @test_warn "more than one random field provided to a geotable, currently we only process the first of these!" check_spatial_data((
        multi_grid,
    ))
end
