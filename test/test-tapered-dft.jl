using SpatialMultitaper, Test
import SpatialMultitaper: tapered_dft, single_tapered_dft, apply_taper, sin_taper_family,
                          DefaultMean, KnownMean
include("test_utilities/TestUtils.jl")
using .TestUtils

@testset "apply_tapers" begin
    @testset "Point data taper application" begin
        # Create test points and marks
        points = PointSet([Point(0.1, 0.2), Point(0.3, 0.4), Point(0.5, 0.6)])
        marks = [1.0, 2.0, 3.0]

        # Create sin taper family
        region = Box(Point(0.0, 0.0), Point(1.0, 1.0))
        taper_family = sin_taper_family((2, 2), region)  # 2x2 = 4 tapers
        tapers = [taper_family[1], taper_family[2]]  # Use first two tapers

        # Apply tapers
        tapered_marks = apply_taper(points, marks, tapers)

        @test length(tapered_marks) == length(points) * length(tapers)
        @test eltype(tapered_marks) == ComplexF64
        @test all(isfinite, tapered_marks)
    end
end

@testset "single_tapered_dft" begin
    using StableRNGs
    rng = StableRNG(42)

    @testset "PointSet (unmarked points)" begin
        # Test with simple point configuration
        points = PointSet([Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0)])
        region = Box(Point(-1.0, -1.0), Point(2.0, 2.0))
        nfreq = (4, 4)
        fmax = (2.0, 2.0)

        # Create sin taper family
        taper_family = sin_taper_family((2, 2), region)  # 2x2 = 4 tapers
        tapers = [taper_family[1], taper_family[2]]  # Use first two tapers

        # Test single_tapered_dft
        result = single_tapered_dft(points, tapers, nfreq, fmax, region, DefaultMean())

        @test size(result) == (4, 4, 2)  # nfreq[1] x nfreq[2] x n_tapers
        @test eltype(result) == ComplexF64
        @test all(isfinite, result)
    end

    @testset "Marked points (GeoTable)" begin
        # Create marked point data
        points, region = make_marked_example(
            rng, n_processes = 1, region_min = (-1.0, -1.0),
            region_max = (2.0, 2.0), point_number = 20)
        data = points[1]

        nfreq = (6, 6)
        fmax = (3.0, 3.0)

        # Create sin taper family
        taper_family = sin_taper_family((2, 2), region)  # 2x2 = 4 tapers
        tapers = [taper_family[1], taper_family[2]]  # Use first two tapers

        result = single_tapered_dft(data, tapers, nfreq, fmax, region, DefaultMean())

        @test size(result) == (6, 6, 2)
        @test eltype(result) == ComplexF64
        @test all(isfinite, result)
    end

    @testset "Gridded data" begin
        # Create gridded test data
        grids, region = make_grids_example(rng, n_processes = 1,
            region_min = (-1.0, -1.0), region_max = (2.0, 2.0),
            grid_dims = (8, 6))
        grid_data = grids[1]

        nfreq = (4, 4)
        # For grid (8,6) with region (-1,-1) to (2,2):
        # spacing = (3/8, 3/6) = (0.375, 0.5)
        # nyquist = (1/(2*0.375), 1/(2*0.5)) = (4/3, 1.0) ≈ (1.333, 1.0)
        # Use multiples: fmax = (4/3, 2.0)
        fmax = (4 / 3, 2.0)

        # Create sin taper family
        taper_family = sin_taper_family((2, 1), region)  # 2x1 = 2 tapers
        tapers = [taper_family[1], taper_family[2]]

        result = single_tapered_dft(domain(grid_data), values(grid_data)[1],
            tapers, nfreq, fmax, region, DefaultMean())

        @test size(result) == (4, 4, 2)
        @test eltype(result) == ComplexF64
        @test all(isfinite, result)
    end
end

@testset "tapered_dft multi-process" begin
    using StableRNGs
    rng = StableRNG(123)

    @testset "Vector input (multiple processes)" begin
        # Create multiple point processes
        points, region = make_points_example(rng, n_processes = 3, as_vector = true,
            region_min = (-1.0, -1.0), region_max = (2.0, 2.0),
            point_number = 25)

        nfreq = (4, 4)
        fmax = (2.0, 2.0)
        # Create sin taper family
        taper_family = sin_taper_family((1, 1), region)  # 1x1 = 1 taper
        tapers = [taper_family[1]]
        mean_method = DefaultMean()  # Single method for all processes

        result = tapered_dft(points, tapers, nfreq, fmax, region, mean_method)

        @test size(result) == (3, 1, 4, 4)  # n_processes x n_tapers x nfreq[1] x nfreq[2]
        @test eltype(result) == ComplexF64
        @test all(isfinite, result)
    end

    @testset "NTuple input (multiple processes)" begin
        # Create data as tuple
        points_data, region = make_points_example(rng, n_processes = 2, as_vector = false,
            region_min = (-1.0, -1.0), region_max = (2.0, 2.0),
            point_number = 30)

        nfreq = (5, 5)
        fmax = (2.5, 2.5)
        # Create sin taper family
        taper_family = sin_taper_family((2, 1), region)  # 2x1 = 2 tapers
        tapers = [taper_family[1], taper_family[2]]
        mean_method = DefaultMean()  # Single method for all processes

        result = tapered_dft(points_data, tapers, nfreq, fmax, region, mean_method)
        @test isa(result, Tuple{Array{ComplexF64, 3}, Array{ComplexF64, 3}})
        @test size(result[1]) == (5, 5, 2)  # nfreq[1] x nfreq[2] x n_tapers
        @test size(result[2]) == (5, 5, 2)
        @test all(isfinite, result[1])
        @test all(isfinite, result[2])
    end

    @testset "Mixed data types" begin
        # Test with mixed point and grid data
        mixed_data, region = make_mixed_example(
            rng, n_processes = (1, 1, 1), as_vector = true,
            region_min = (-0.5, -0.5), region_max = (1.5, 1.5),
            grid_dims = (6, 6), point_number = 20)

        nfreq = (4, 4)
        # For grid (6,6) with region (-0.5,-0.5) to (1.5,1.5):
        # spacing = (2/6, 2/6) = (1/3, 1/3) ≈ (0.333, 0.333)
        # nyquist = (1/(2*0.333), 1/(2*0.333)) = (1.5, 1.5)
        # Use multiples: fmax = (1.5, 3.0) = (1, 2) * nyquist
        fmax = (1.5, 3.0)
        # Create sin taper family
        taper_family = sin_taper_family((1, 1), region)  # 1x1 = 1 taper
        tapers = [taper_family[1]]
        mean_method = DefaultMean()  # Single method for all processes

        result = tapered_dft(mixed_data, tapers, nfreq, fmax, region, mean_method)

        @test size(result) == (3, 1, 4, 4)  # 3 processes (point, grid, marked) x 1 taper x 4 x 4 freqs
        @test eltype(result) == ComplexF64
        @test all(isfinite, result)
    end
end

@testset "Different mean methods" begin
    using StableRNGs
    rng = StableRNG(456)

    # Create test data
    marked_data, region = make_marked_example(rng, n_processes = 1,
        region_min = (-1.0, -1.0), region_max = (1.0, 1.0),
        point_number = 50)
    data = marked_data[1]

    nfreq = (6, 6)
    fmax = (3.0, 3.0)
    # Create sin taper family
    taper_family = sin_taper_family((1, 1), region)  # 1x1 = 1 taper
    tapers = [taper_family[1]]

    @testset "DefaultMean vs KnownMean" begin
        result_default_mean = single_tapered_dft(
            data, tapers, nfreq, fmax, region, DefaultMean())
        result_known_mean = single_tapered_dft(
            data, tapers, nfreq, fmax, region, KnownMean(0.5))

        @test size(result_default_mean) ==
              size(result_known_mean) == (6, 6, 1)
        @test eltype(result_default_mean) ==
              eltype(result_known_mean) == ComplexF64

        # Results should be different across mean methods
        @test !all(result_default_mean .≈ result_known_mean)

        # All should be finite
        @test all(isfinite, result_default_mean)
        @test all(isfinite, result_known_mean)
    end
end

@testset "Edge cases and error handling" begin
    using StableRNGs
    rng = StableRNG(789)

    # @testset "Empty data" begin # TODO: need to catch this
    #     points = PointSet(typeof(Point(1, 1))[])
    #     region = Box(Point(-1.0, -1.0), Point(1.0, 1.0))
    #     nfreq = (4, 4)
    #     fmax = (2.0, 2.0)
    #     # Create sin taper family
    #     taper_family = sin_taper_family((1, 1), region)  # 1x1 = 1 taper
    #     tapers = [taper_family[1]]

    #     result = single_tapered_dft(points, tapers, nfreq, fmax, region, DefaultMean())
    #     @test size(result) == (4, 4, 1)
    #     @test eltype(result) == ComplexF64
    # end

    @testset "Single point" begin
        points = PointSet([Point(0.0, 0.0)])
        region = Box(Point(-1.0, -1.0), Point(1.0, 1.0))
        nfreq = (3, 3)
        fmax = (1.5, 1.5)
        # Create sin taper family
        taper_family = sin_taper_family((1, 1), region)  # 1x1 = 1 taper
        tapers = [taper_family[1]]

        result = single_tapered_dft(points, tapers, nfreq, fmax, region, DefaultMean())
        @test size(result) == (3, 3, 1)
        @test all(isfinite, result)
    end

    @testset "Multiple tapers" begin
        points, region = make_points_example(rng, n_processes = 1,
            region_min = (-1.0, -1.0), region_max = (1.0, 1.0))
        data = points[1]

        nfreq = (4, 4)
        fmax = (2.0, 2.0)

        # Create sin taper family with multiple tapers
        taper_family = sin_taper_family((3, 1), region)  # 3x1 = 3 tapers
        tapers = [taper_family[1], taper_family[2], taper_family[3]]

        result = single_tapered_dft(data, tapers, nfreq, fmax, region, DefaultMean())
        @test size(result) == (4, 4, 3)  # 3 tapers
        @test all(isfinite, result)
    end
end

@testset "1D cases" begin
    using StableRNGs
    rng = StableRNG(321)

    @testset "1D point process" begin
        points, region = make_points_example(rng, n_processes = 1, as_vector = true,
            region_min = (-2.0,), region_max = (3.0,), point_number = 20)
        data = points[1]

        nfreq = (8,)
        fmax = (4.0,)
        # Create sin taper family for 1D
        taper_family = sin_taper_family((2,), region)  # 2 tapers in 1D
        tapers = [taper_family[1]]

        result = single_tapered_dft(data, tapers, nfreq, fmax, region, DefaultMean())
        @test size(result) == (8, 1)  # 1D frequency space x 1 taper
        @test all(isfinite, result)
    end

    @testset "1D grid" begin
        grids, region = make_grids_example(rng, n_processes = 1, as_vector = true,
            region_min = (-1.0,), region_max = (2.0,), grid_dims = (12,))
        grid_data = grids[1]

        nfreq = (6,)
        # For grid (12,) with region (-1,) to (2,):
        # spacing = 3/12 = 0.25
        # nyquist = 1/(2*0.25) = 2.0
        # Use multiple: fmax = (4.0,) = 2 * nyquist
        fmax = (4.0,)
        # Create sin taper family for 1D
        taper_family = sin_taper_family((2,), region)  # 2 tapers in 1D
        tapers = [taper_family[1]]

        result = single_tapered_dft(domain(grid_data), values(grid_data)[1],
            tapers, nfreq, fmax, region, DefaultMean())
        @test size(result) == (6, 1)
        @test all(isfinite, result)
    end
end
