using SpatialMultitaper, Test
import SpatialMultitaper: intensity_index, PoissonProcess, PartialResampler, StandardShift,
                          UniformShift, getestimate
import Random

@testset "ToroidalShift" begin
    region = Box(Point(0, 0), Point(100, 100))
    shift = ToroidalShift(region)
    @test all(shift.shift.min .≈ (-50.0, -50.0))
    @test all(shift.shift.max .≈ (50.0, 50.0))
    shift_0 = rand(shift.shift)
    @test all(x -> -50 ≤ x ≤ 50, shift_0)
end

@testset "intensity_index" begin
    point = Point(5.23)
    grid = CartesianGrid((10,), Point(5.0), (0.1,))
    @test intensity_index(point, minimum(grid), spacing(grid)) == CartesianIndex(3)
end

@testset "partial K resampling" begin
    region = Box(Point(0, 0), Point(100, 100))
    pattern = rand(PoissonProcess(0.01), region)
    pattern2 = rand(PoissonProcess(0.01), region)
    pattern3 = rand(PoissonProcess(0.01), region)
    tapers = sin_taper_family((4, 4), region)
    nfreq = (100, 100)
    fmax = (nfreq .- 1) ./ (2 .* 100)
    radii = 0.3:0.1:0.5
    for data in [(pattern, pattern2), (pattern, pattern2, pattern3)]
        resampler = PartialResampler(
            data,
            region;
            tapers = tapers,
            nfreq = nfreq,
            fmax = fmax,
            nfreq_marginal_compute = (50, 50),
            fmax_marginal_compute = (0.2, 0.2),
            shift_method = ToroidalShift(region)
        )

        results = partial_shift_resample(
            partial_k_function,
            resampler;
            radii = radii
        )
        @test results.radii == radii

        resampler_2 = PartialResampler(
            data,
            region;
            tapers = tapers,
            nfreq = nfreq,
            fmax = fmax,
            nfreq_marginal_compute = (50, 50),
            fmax_marginal_compute = (0.2, 0.2),
            shift_method = StandardShift(UniformShift((-0.1, -0.1), (0.1, 0.1)))
        )

        results = partial_shift_resample(
            partial_k_function,
            resampler_2;
            radii = radii
        )
        @test results.radii == radii

        rng = Random.MersenneTwister(1234)
        resampled_data_1 = rand(rng, resampler.marginal_resampler[1])
        rng = Random.MersenneTwister(1234)
        resampled_data_2 = rand(rng, resampler.marginal_resampler[1])
        @test observations(resampled_data_1) == observations(resampled_data_2)

        rng = Random.MersenneTwister(1234)
        results1 = shift_resample(
            rng,
            spatial_data(data, region),
            k_function,
            ToroidalShift(region);
            radii = radii,
            tapers = tapers,
            nfreq = nfreq,
            fmax = fmax
        )

        rng = Random.MersenneTwister(1234)
        results2 = shift_resample(
            rng,
            spatial_data(data, region),
            k_function,
            ToroidalShift(region);
            radii = radii,
            tapers = tapers,
            nfreq = nfreq,
            fmax = fmax
        )
        results3 = shift_resample(
            rng,
            spatial_data(data, region),
            k_function,
            ToroidalShift(region);
            radii = radii,
            tapers = tapers,
            nfreq = nfreq,
            fmax = fmax
        )
        @test getestimate(results1) == getestimate(results2)
        @test getestimate(results2) !== getestimate(results3)

        rng = Random.MersenneTwister(1234)
        results1 = partial_shift_resample(
            rng,
            partial_k_function,
            resampler;
            radii = radii
        )
        rng = Random.MersenneTwister(1234)
        results2 = partial_shift_resample(
            rng,
            partial_k_function,
            resampler;
            radii = radii
        )
        results3 = partial_shift_resample(
            rng,
            partial_k_function,
            resampler;
            radii = radii
        )
        @test getestimate(results1) == getestimate(results2)
        @test getestimate(results2) !== getestimate(results3)
    end
end
