@testset "K function" begin
    region = Box(Point(0,0), Point(3,3))
    pattern = PointSet([Point(0,0), Point(1,1)])
    pattern2 = PointSet([Point(0.3,0.2), Point(0.8,0.4), Point(0.5,0.5)])
    griddata = georef((rf=rand(6*6),), CartesianGrid((0,0), (3,3), dims=(6,6)))
    tapers = sin_taper_family((3,3), region)
    nfreq = (10,10)
    fmax = (2,2)
    data = (pattern, pattern2, griddata)
    R = 0.3:0.1:0.5
    results = partial_K(data, R, tapers; region=region, nfreq=nfreq, fmax=fmax)
    @test results[1] == R
    @test results[2] isa Dict
end