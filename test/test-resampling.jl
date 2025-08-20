@testset "ToroidalShift" begin
    region = Box(Point(0, 0), Point(100, 100))
    shift = ToroidalShift(region)
    @test all(shift.shift.min .≈ (-50.0, -50.0))
    @test all(shift.shift.max .≈ (50.0, 50.0))
    shift_0 = rand(shift.shift)
    @test all(x -> -50 ≤ x ≤ 50, shift_0)
end

@testset "partial K resampling" begin
    region = Box(Point(0, 0), Point(3, 3))
    pattern = PointSet([Point(0, 0), Point(1, 1)])
    pattern2 = PointSet([Point(0.3, 0.2), Point(0.8, 0.4), Point(0.5, 0.5)])
    pattern3 = PointSet([Point(0.1, 0.1), Point(0.2, 0.2), Point(0.3, 0.3)])
    tapers = sin_taper_family((3, 3), region)
    nfreq = (10, 10)
    fmax = (2, 2)
    data = (pattern, pattern2, pattern3)
    radii = 0.3:0.1:0.5
    results = partial_shift_resample(
        data,
        region,
        partial_K_function,
        ToroidalShift(region);
        radii = radii,
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
    )
    @test results.radii == radii

    results = partial_shift_resample(
        data,
        region,
        partial_K_function,
        Spmt.StandardShift(Spmt.UniformShift((-0.1, -0.1), (0.1, 0.1)));
        radii = radii,
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
    )
    @test results.radii == radii
end
