using SpatialMultitaper, Test, StableRNGs
include("../test_utilities/TestUtils.jl")
using .TestUtils

import SpatialMultitaper: CFunction, KFunction, c_function, k_function, is_partial,
                          partial_c_function, partial_k_function, sdf2C, C2K, getestimate,
                          getargument, getbaseestimatename, ProcessInformation,
                          getestimationinformation, getprocessinformation,
                          EstimationInformation, MarginalTrait,
                          mean_estimate, DefaultMean

@testset "C Function" begin
    rng = StableRNG(123)

    @testset "CFunction struct" begin
        radii = [0.1, 0.2, 0.3]
        values = [1.0, 0.8, 0.6]
        processinfo = ProcessInformation([1], [1], ones(1, 1), ones(1, 1), Val{2}())
        estimationinfo = EstimationInformation(5)

        c_func = CFunction{MarginalTrait}(radii, values, processinfo, estimationinfo)

        @test getargument(c_func) == radii
        @test getestimate(c_func) == values
        @test embeddim(c_func) == 2
        @test size(c_func) == (1, 1)
        @test getbaseestimatename(CFunction) == "C function"
    end

    @testset "c_function from data" begin
        data, region = make_points_example(rng, n_processes = 2, point_number = 50)
        radii = [0.05, 0.1, 0.15, 0.2]

        c_est = c_function(data, region,
            radii = radii,
            nfreq = (8, 8),
            fmax = (0.4, 0.4),
            tapers = sin_taper_family((2, 2), region))

        @test c_est isa CFunction
        @test getargument(c_est) == radii
        @test length(getestimate(c_est)) == length(radii)
        @test size(c_est) == (2, 2)
    end

    @testset "c_function from spectra" begin
        data, region = make_points_example(rng, n_processes = 2, point_number = 40)
        spec = spectra(data, region, nfreq = (6, 6), fmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        radii = [0.1, 0.2, 0.3]
        c_est = c_function(spec, radii = radii)

        @test c_est isa CFunction
        @test getargument(c_est) == radii
        @test getprocessinformation(c_est) == getprocessinformation(spec)
        @test getestimationinformation(c_est) == getestimationinformation(spec)
    end

    @testset "partial_c_function" begin
        data, region = make_points_example(rng, n_processes = 3, point_number = 40)
        radii = [0.05, 0.1, 0.15]

        c_partial = partial_c_function(data, region,
            radii = radii,
            nfreq = (6, 6),
            fmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        @test c_partial isa CFunction
        @test is_partial(c_partial) == true
        @test getargument(c_partial) == radii
    end
end

@testset "K Function" begin
    rng = StableRNG(123)

    @testset "KFunction struct" begin
        radii = [0.1, 0.2, 0.3]
        values = [0.01, 0.04, 0.09]  # Roughly π*r² for circle
        processinfo = ProcessInformation([1], [1], ones(1, 1), ones(1, 1), Val{2}())
        estimationinfo = EstimationInformation(5)

        k_func = KFunction{MarginalTrait}(radii, values, processinfo, estimationinfo)

        @test getargument(k_func) == radii
        @test getestimate(k_func) == values
        @test getbaseestimatename(KFunction) == "K function"
    end

    @testset "k_function from c_function" begin
        data, region = make_points_example(rng, n_processes = 2, point_number = 40)
        radii = [0.05, 0.1, 0.15]

        c_est = c_function(data, region,
            radii = radii,
            nfreq = (6, 6),
            fmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        k_est = k_function(c_est)

        @test k_est isa KFunction
        @test getargument(k_est) == getargument(c_est)
        @test size(k_est) == size(c_est)
        @test getprocessinformation(k_est) == getprocessinformation(c_est)
    end

    @testset "k_function from data" begin
        data, region = make_points_example(rng, n_processes = 2, point_number = 40)
        radii = [0.05, 0.1, 0.15]

        k_est = k_function(data, region,
            radii = radii,
            nfreq = (6, 6),
            fmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        @test k_est isa KFunction
        @test getargument(k_est) == radii
    end
end

@testset "Mathematical Properties" begin
    rng = StableRNG(123)

    @testset "C function properties" begin
        # C function should be non-negative for typical point processes
        data, region = make_points_example(rng, n_processes = 1, point_number = 50)
        radii = [0.05, 0.1, 0.15, 0.2]

        c_est = c_function(data, region,
            radii = radii,
            nfreq = (8, 8),
            fmax = (0.4, 0.4),
            tapers = sin_taper_family((2, 2), region))

        values = getestimate(c_est)
        # For typical point processes, C function is often positive
        # but this depends on the specific process and may have negative values
        @test all(isfinite, values)
        @test length(values) == length(radii)
    end

    @testset "K function properties" begin
        # K function should generally be monotonic increasing for clustered processes
        rng = StableRNG(123)
        data, region = make_points_example(rng, n_processes = 1, point_number = 60)
        radii = [0.05, 0.1, 0.15, 0.2, 0.25]

        k_est = k_function(data, region,
            radii = radii,
            nfreq = (10, 10),
            fmax = (0.5, 0.5),
            tapers = sin_taper_family((2, 2), region))

        values = getestimate(k_est)
        @test all(isfinite, values)
        @test all(values .≥ 0)  # K function should be non-negative

        # For complete spatial randomness, K(r) ≈ π*r² in 2D
        # Test that values are reasonable order of magnitude
        expected_orders = π .* radii .^ 2

        @test all(values ./ expected_orders .> 0.1)  # Within order of magnitude
        @test all(values ./ expected_orders .< 10.0)
    end
end

@testset "Indexing and Access" begin
    rng = StableRNG(123)
    data, region = make_points_example(rng, n_processes = 2, point_number = 30)
    radii = [0.1, 0.2, 0.3]

    c_est = c_function(data, region,
        radii = radii,
        nfreq = (4, 4),
        fmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    @testset "Process indexing" begin
        c_sub = c_est[1, 2]
        @test c_sub isa CFunction
        @test size(c_sub) == (1, 1)
        @test getargument(c_sub) == getargument(c_est)
    end

    @testset "Radius indexing" begin
        c_radius = c_est[1, 1, 2]  # Second radius
        @test c_radius isa CFunction
        @test length(getargument(c_radius)) == 1
        @test getargument(c_radius)[1] ≈ radii[2]
    end
end

@testset "Integration with Other Estimates" begin
    rng = StableRNG(123)
    data, region = make_points_example(rng, n_processes = 2, point_number = 40)

    @testset "Partial spatial functions" begin
        radii = [0.05, 0.1, 0.15]

        # Test partial C function
        c_partial = partial_c_function(data, region,
            radii = radii,
            nfreq = (6, 6),
            fmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        @test is_partial(c_partial) == true

        # Test partial K function
        k_partial = partial_k_function(data, region,
            radii = radii,
            nfreq = (6, 6),
            fmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        @test is_partial(k_partial) == true
    end
end

@testset "Edge Cases" begin
    rng = StableRNG(123)

    @testset "Single radius" begin
        data, region = make_points_example(rng, n_processes = 1, point_number = 30)

        c_est = c_function(data, region,
            radii = [0.1],
            nfreq = (4, 4),
            fmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        @test length(getargument(c_est)) == 1
        @test length(getestimate(c_est)) == 1
    end

    @testset "Very small radii" begin
        data, region = make_points_example(rng, n_processes = 1, point_number = 30)

        c_est = c_function(data, region,
            radii = [0.001, 0.01],
            nfreq = (6, 6),
            fmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        values = getestimate(c_est)
        @test all(isfinite, values)
    end
end
