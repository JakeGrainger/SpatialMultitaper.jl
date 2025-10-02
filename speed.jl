##
using SpatialMultitaper
function foo()
    points1 = PointSet(Point.(rand(1000), rand(1000)))
    points2 = PointSet(Point.(rand(1000), rand(1000)))
    points3 = PointSet(Point.(rand(1000), rand(1000)))
    box = Box((0, 0), (1, 1))
    data = spatial_data((points1, points2, points3), box)
    radii = range(0, 0.25, 100)
    tapers = sin_taper_family((4, 4), box)
    @profview partial_l_function(data, radii = radii, kmax = 100, tapers = tapers)
end

foo()
