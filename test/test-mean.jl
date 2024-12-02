region = Box(Point(0, 0), Point(2, 2))
points = PointSet([Point(0, 0), Point(2, 2), Point(0.5, 0.5), Point(2.5, 0.5)])
marks = [1, 2, 3, 4]
rf = [1, 2, 3, 4, 5, 6, 7, 8, 9]
grid = CartesianGrid((0.0, 0.0), (3.0, 3.0), dims = (3, 3))
pointdata = georef((marks = marks,), points)
griddata = georef((rf = rf,), grid)

@test SpatialMultitaper.mean_estimate(pointdata, region, DefaultMean()) == 3 / 2 # (1+2+3)/4 (because of the region, and that the region has measure 4)
@test SpatialMultitaper.mean_estimate(griddata, region, DefaultMean()) == 3.0 # (1+2+4+5)/4 (because of the region)
@test SpatialMultitaper.mean_estimate(pointdata, points, region, KnownMean(1.0)) == 1.0
@test SpatialMultitaper.mean_estimate(griddata, grid, region, KnownMean(1.0)) == 1.0
