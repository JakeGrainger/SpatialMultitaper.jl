@testset "ToroidalShift" begin
    region = Box(Point(0, 0), Point(100, 100))
    shift = ToroidalShift(region)
    @test all(shift.shift.min .≈ (-50.0, -50.0))
    @test all(shift.shift.max .≈ (50.0, 50.0))
    shift_0 = rand(shift.shift)
    @test all(x -> -50 ≤ x ≤ 50, shift_0)
end

@testset "partial K resampling" begin
    region = Box(Point(0, 0), Point(100, 100))
    pattern = rand(Spmt.PoissonProcess(0.01), region)
    pattern2 = rand(Spmt.PoissonProcess(0.01), region)
    pattern3 = rand(Spmt.PoissonProcess(0.01), region)
    tapers = sin_taper_family((4, 4), region)
    nfreq = (100, 100)
    fmax = (2, 2)
    data = (pattern, pattern2, pattern3)
    radii = 0.3:0.1:0.5

    resampler = Spmt.MarginalResampler(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        radii = 0:0.01:2,
        grid = CartesianGrid(Point(-1, -1), Point(4, 4), dims = (100, 100)),
    )

    results = partial_shift_resample(
        data,
        region,
        partial_K_function,
        ToroidalShift(region),
        resampler;
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
        Spmt.StandardShift(Spmt.UniformShift((-0.1, -0.1), (0.1, 0.1))),
        resampler;
        radii = radii,
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
    )
    @test results.radii == radii

    rng = Spmt.Random.MersenneTwister(1234)
    resampled_data_1 = rand(rng, resampler[1])
    rng = Spmt.Random.MersenneTwister(1234)
    resampled_data_2 = rand(rng, resampler[1])
    @test resampled_data_1 == resampled_data_2

    rng = Spmt.Random.MersenneTwister(1234)
    results1 = shift_resample(
        rng,
        data,
        region,
        K_function,
        ToroidalShift(region);
        radii = radii,
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
    )

    rng = Spmt.Random.MersenneTwister(1234)
    results2 = shift_resample(
        rng,
        data,
        region,
        K_function,
        ToroidalShift(region);
        radii = radii,
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
    )
    @test results1.K_function == results2.K_function

    rng = Spmt.Random.MersenneTwister(1234)
    results1 = partial_shift_resample(
        rng,
        data,
        region,
        partial_K_function,
        ToroidalShift(region),
        resampler;
        radii = radii,
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
    )
    rng = Spmt.Random.MersenneTwister(1234)
    results2 = partial_shift_resample(
        rng,
        data,
        region,
        partial_K_function,
        ToroidalShift(region),
        resampler;
        radii = radii,
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
    )
    @test results1.partial_K_function == results2.partial_K_function
end
