import SpatialMultitaper: Box, Point
@testset "partial_covariance_density" begin
    region = Box(Point(0,0), Point(3,3))
    pattern = georef((marks=[1,2.4],), PointSet([Point(0,0), Point(1,1)]))
    griddata = georef((rf=rand(6*6),), CartesianGrid((0,0), (3,3), dims=(6,6)))
    tapers = sin_taper_family((3,3), region)
    nfreq = (10,10)
    fmax = (2,2)
    mt_est = multitaper_estimate((pattern,griddata), tapers, nfreq=nfreq, fmax=fmax, region=region)
    @test partial_covariance_density(mt_est) isa NamedTuple
end