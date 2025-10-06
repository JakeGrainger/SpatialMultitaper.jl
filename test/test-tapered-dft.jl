using SpatialMultitaper, Test, StaticArrays
import SpatialMultitaper: tapered_dft, _single_tapered_dft, _apply_taper!,
                          DefaultMean, KnownMean
include("test_utilities/TestData.jl")
using .TestData

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
        tapered_marks = zeros(ComplexF64, length(points) * length(tapers))
        _apply_taper!(tapered_marks, points, marks, tapers)

        @test length(tapered_marks) == length(points) * length(tapers)
        @test eltype(tapered_marks) == ComplexF64
        @test all(isfinite, tapered_marks)
    end
end

@testset "_single_tapered_dft" begin
    using StableRNGs
    rng = StableRNG(42)

    @testset "PointSet (unmarked points)" begin
        # Test with simple point configuration
        points = PointSet([Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0)])
        region = Box(Point(-1.0, -1.0), Point(2.0, 2.0))
        nk = (4, 4)
        kmax = (2.0, 2.0)

        # Create sin taper family
        taper_family = sin_taper_family((2, 2), region)  # 2x2 = 4 tapers
        tapers = [taper_family[1], taper_family[2]]  # Use first two tapers

        # Test _single_tapered_dft
        result = _single_tapered_dft(
            spatial_data(points, region), tapers, nk, kmax, DefaultMean())

        @test size(result) == (4, 4, 2)  # nk[1] x nk[2] x n_tapers
        @test eltype(result) == ComplexF64
        @test all(isfinite, result)
    end

    @testset "Marked points (GeoTable)" begin
        # Create marked point data
        marked_data = make_marked_example(
            rng, n_processes = 1, return_type = :single,
            region_min = (-1.0, -1.0), region_max = (2.0, 2.0), point_number = 20)

        nk = (6, 6)
        kmax = (3.0, 3.0)

        # Create sin taper family
        region = getregion(marked_data)
        tapers = sin_taper_family((2, 2), region)  # 2x2 = 4 tapers

        result = _single_tapered_dft(marked_data, tapers, nk, kmax, DefaultMean())

        @test size(result) == (6, 6, 4)
        @test eltype(result) == ComplexF64
        @test all(isfinite, result)
    end

    @testset "Gridded data" begin
        # Create gridded test data
        grid_data = make_grids_example(rng, n_processes = 1, return_type = :single,
            region_min = (-1.0, -1.0), region_max = (2.0, 2.0), grid_dims = (8, 6))

        nk = (4, 4)
        # For grid (8,6) with region (-1,-1) to (2,2):
        # spacing = (3/8, 3/6) = (0.375, 0.5)
        # nyquist = (1/(2*0.375), 1/(2*0.5)) = (4/3, 1.0) ≈ (1.333, 1.0)
        # Use multiples: kmax = (4/3, 2.0)
        kmax = (4 / 3, 2.0)

        # Create sin taper family
        region = getregion(grid_data)
        tapers = sin_taper_family((2, 1), region)  # 2x1 = 2 tapers

        result = _single_tapered_dft(grid_data, tapers, nk, kmax, DefaultMean())

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
        points_data = make_points_example(rng, n_processes = 3, return_type = :vector,
            region_min = (-1.0, -1.0), region_max = (2.0, 2.0), point_number = 25)

        nk = (4, 4)
        kmax = (2.0, 2.0)
        # Create sin taper family
        region = getregion(points_data)
        tapers = sin_taper_family((1, 1), region)  # 1x1 = 1 taper
        mean_method = DefaultMean()  # Single method for all processes

        result = tapered_dft(points_data, tapers, nk, kmax, mean_method)

        @test size(result) == (3, 1, 4, 4)  # n_processes x n_tapers x nk[1] x nk[2]
        @test eltype(result) == ComplexF64
        @test all(isfinite, result)
    end

    @testset "NTuple input (multiple processes)" begin
        # Create data as tuple
        points_data = make_points_example(rng, n_processes = 2, return_type = :tuple,
            region_min = (-1.0, -1.0), region_max = (2.0, 2.0), point_number = 30)

        nk = (5, 5)
        kmax = (2.5, 2.5)
        # Create sin taper family
        region = getregion(points_data)
        tapers = sin_taper_family((3, 1), region)
        mean_method = DefaultMean()  # Single method for all processes

        result = tapered_dft(points_data, tapers, nk, kmax, mean_method)
        @test isa(result, Array{SVector{2, ComplexF64}, 3}) # 2 processes, 2 dims in space and then tapers
        @test size(result) == (5, 5, 3)  # nk[1] x nk[2] x n_tapers

        @test all(x -> all(isfinite.(x)), result)
    end

    @testset "Mixed data types" begin
        # Test with mixed point and grid data
        mixed_data = make_mixed_example(
            rng, n_processes = (1, 1, 1), return_type = :vector,
            region_min = (-0.5, -0.5), region_max = (1.5, 1.5),
            grid_dims = (6, 6), point_number = 20)

        nk = (4, 4)
        # For grid (6,6) with region (-0.5,-0.5) to (1.5,1.5):
        # spacing = (2/6, 2/6) = (1/3, 1/3) ≈ (0.333, 0.333)
        # nyquist = (1/(2*0.333), 1/(2*0.333)) = (1.5, 1.5)
        # Use multiples: kmax = (1.5, 3.0) = (1, 2) * nyquist
        kmax = (1.5, 3.0)
        # Create sin taper family
        region = getregion(mixed_data)
        tapers = sin_taper_family((1, 1), region)  # 1x1 = 1 taper
        mean_method = DefaultMean()  # Single method for all processes

        result = tapered_dft(mixed_data, tapers, nk, kmax, mean_method)

        @test size(result) == (3, 1, 4, 4)  # 3 processes (point, grid, marked) x 1 taper x 4 x 4 ks
        @test eltype(result) == ComplexF64
        @test all(isfinite, result)
    end
end

@testset "Different mean methods" begin
    using StableRNGs
    rng = StableRNG(456)

    # Create test data
    marked_data = make_marked_example(rng, n_processes = 1, return_type = :single,
        region_min = (-1.0, -1.0), region_max = (1.0, 1.0), point_number = 50)

    nk = (6, 6)
    kmax = (3.0, 3.0)
    # Create sin taper family
    region = getregion(marked_data)
    tapers = sin_taper_family((1, 1), region)  # 1x1 = 1 taper

    @testset "DefaultMean vs KnownMean" begin
        result_default_mean = _single_tapered_dft(
            marked_data, tapers, nk, kmax, DefaultMean())
        result_known_mean = _single_tapered_dft(
            marked_data, tapers, nk, kmax, KnownMean(0.5))

        @test size(result_default_mean) == size(result_known_mean) == (6, 6, 1)
        @test eltype(result_default_mean) == eltype(result_known_mean) == ComplexF64

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
    #     nk = (4, 4)
    #     kmax = (2.0, 2.0)
    #     # Create sin taper family
    #     taper_family = sin_taper_family((1, 1), region)  # 1x1 = 1 taper
    #     tapers = [taper_family[1]]

    #     result = _single_tapered_dft(points, tapers, nk, kmax, region, DefaultMean())
    #     @test size(result) == (4, 4, 1)
    #     @test eltype(result) == ComplexF64
    # end

    @testset "Single point" begin
        points = PointSet([Point(0.0, 0.0)])
        region = Box(Point(-1.0, -1.0), Point(1.0, 1.0))
        nk = (3, 3)
        kmax = (1.5, 1.5)
        # Create sin taper family
        tapers = sin_taper_family((1, 1), region)  # 1x1 = 1 taper

        result = _single_tapered_dft(
            spatial_data(points, region), tapers, nk, kmax, DefaultMean())
        @test size(result) == (3, 3, 1)
        @test all(isfinite, result)
    end

    @testset "Multiple tapers" begin
        points_data = make_points_example(rng, n_processes = 1, return_type = :single,
            region_min = (-1.0, -1.0), region_max = (1.0, 1.0))

        nk = (4, 4)
        kmax = (2.0, 2.0)

        # Create sin taper family with multiple tapers
        region = getregion(points_data)
        tapers = sin_taper_family((3, 1), region)  # 3x1 = 3 tapers

        result = _single_tapered_dft(points_data, tapers, nk, kmax, DefaultMean())
        @test size(result) == (4, 4, 3)  # 3 tapers
        @test all(isfinite, result)
    end
end

@testset "1D cases" begin
    using StableRNGs
    rng = StableRNG(321)

    @testset "1D point process" begin
        points_data = make_points_example(rng, n_processes = 1, return_type = :single,
            dim = 1, region_min = (-2.0,), region_max = (3.0,), point_number = 20)

        nk = (8,)
        kmax = (4.0,)
        # Create sin taper family for 1D
        region = getregion(points_data)
        tapers = sin_taper_family((1,), region)  # 1 taper in 1D

        result = _single_tapered_dft(points_data, tapers, nk, kmax, DefaultMean())
        @test size(result) == (8, 1)  # 1D wavenumber space x 1 taper
        @test all(isfinite, result)
    end

    @testset "1D grid" begin
        grid_data = make_grids_example(rng, n_processes = 1, return_type = :single,
            dim = 1, region_min = (-1.0,), region_max = (2.0,), grid_dims = (12,))

        nk = (6,)
        # For grid (12,) with region (-1,) to (2,):
        # spacing = 3/12 = 0.25
        # nyquist = 1/(2*0.25) = 2.0
        # Use multiple: kmax = (4.0,) = 2 * nyquist
        kmax = (4.0,)
        # Create sin taper family for 1D
        region = getregion(grid_data)
        tapers = sin_taper_family((2,), region)  # 2 tapers in 1D

        result = _single_tapered_dft(grid_data, tapers, nk, kmax, DefaultMean())
        @test size(result) == (6, 2)
        @test all(isfinite, result)
    end
end
