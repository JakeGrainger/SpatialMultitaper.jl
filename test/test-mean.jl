using SpatialMultitaper, Test, StaticArrays

import SpatialMultitaper: mean_estimate

region = Box(Point(0, 0), Point(2, 2))

points = PointSet([Point(0, 0), Point(2, 2), Point(0.5, 0.5), Point(2.5, 0.5)])
marks = [1, 2, 3, 4]
pointdata = georef((marks = marks,), points)
@test mean_estimate(pointdata, region, DefaultMean()) == 3 / 2 # (1+2+3)/4 (because of the region, and that the region has measure 4)
@test mean_estimate(points, region, DefaultMean()) == 3 / 4 # because the region has measure 4 and there are 3 points in the region
@test mean_estimate(pointdata, region, DefaultMean()) isa Number

grid = CartesianGrid((0.0, 0.0), (3.0, 3.0), dims = (3, 3))
rf = [1, 2, 3, 4, 5, 6, 7, 8, 9]
griddata = georef((rf = rf,), grid)
@test mean_estimate(griddata, region, DefaultMean()) == 3.0 # (1+2+4+5)/4 (because not all the points are in the region)

@test mean_estimate(pointdata, region, KnownMean(1.0)) == 1.0
@test mean_estimate(griddata, region, KnownMean(1.0)) == 1.0

@test mean_estimate(spatial_data((points,), region), DefaultMean()) isa SVector{1}
