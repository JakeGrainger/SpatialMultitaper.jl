using Test, SpatialMultitaper, StaticArrays, LinearAlgebra, StableRNGs

import SpatialMultitaper: covariance_zero_atom, unitless_measure

# Import test utilities
include("test_utilities/TestData.jl")
using .TestData

@testset "covariance_zero_atom tests" begin
    rng = StableRNG(123)

    @testset "PointSet tests" begin
        points_data, region = make_points_example(rng, n_processes = 1, point_number = 100)

        # Test basic computation
        result = covariance_zero_atom(spatial_data(points_data, region))
        expected = length(points_data[1]) / unitless_measure(region)
        @test result[1] ≈ expected
        @test result[1] > 0  # Should be positive for non-empty point set

        # Test with different point densities
        sparse_data, _ = make_points_example(rng, n_processes = 1, point_number = 10)
        sparse_result = covariance_zero_atom(spatial_data(sparse_data[1], region))
        @test sparse_result < result[1]  # Lower density should give smaller atom
    end

    @testset "Marked point tests (GeoTable)" begin
        marked_data, region = make_marked_example(rng, n_processes = 1, point_number = 50)

        result = covariance_zero_atom(spatial_data(marked_data, region))

        # Should compute sum of squared marks divided by region measure
        marks = values(marked_data[1])[1]  # Get the marks
        expected = sum(abs2, marks) / unitless_measure(region)
        @test result[1] ≈ expected
        @test result[1] > 0  # Should be positive for non-zero marks
    end

    @testset "CartesianGrid tests" begin
        grids_data, region = make_grids_example(rng, n_processes = 1, grid_dims = (10, 8))

        # Grid data should always return zero atom
        result = covariance_zero_atom(spatial_data(grids_data, region))
        @test result == SMatrix{1, 1}(0.0)

        # Test with different grid sizes
        large_grids, _ = make_grids_example(rng, n_processes = 1, grid_dims = (50, 30))
        large_result = covariance_zero_atom(spatial_data(large_grids, region))
        @test large_result == SMatrix{1, 1}(0.0)
    end

    @testset "Multiple processes - Tuple format" begin
        # Test with multiple point processes
        points_data, region = make_points_example(rng, n_processes = 3, point_number = 60)

        result = covariance_zero_atom(spatial_data(points_data, region))
        @test result isa AbstractMatrix
        @test size(result) == (3, 3)
        @test isdiag(result)  # Should be diagonal matrix

        # Diagonal elements should be individual atom estimates
        for i in 1:3
            individual_result = covariance_zero_atom(spatial_data(points_data[i], region))
            @test result[i, i] ≈ individual_result
            @test result[i, i] > 0
        end
    end

    @testset "Multiple processes - Vector format" begin
        points_data, region = make_points_example(rng, n_processes = 2, as_vector = true)

        result = covariance_zero_atom(spatial_data(points_data, region))
        @test result isa AbstractMatrix
        @test size(result) == (2, 2)
        @test isdiag(result)

        # Compare with individual computations
        individual_1 = covariance_zero_atom(spatial_data(points_data[1], region))
        individual_2 = covariance_zero_atom(spatial_data(points_data[2], region))
        @test result[1, 1] ≈ individual_1
        @test result[2, 2] ≈ individual_2
    end

    @testset "Single element tuple" begin
        points_data, region = make_points_example(rng, n_processes = 1)

        result = covariance_zero_atom(spatial_data(points_data, region))
        expected = covariance_zero_atom(spatial_data(points_data[1], region))
        @test result[1] ≈ expected
        @test result isa SMatrix
        @test expected isa Number
    end

    @testset "Mixed data types" begin
        mixed_data, region = make_mixed_example(rng, n_processes = (1, 1, 1))

        result = covariance_zero_atom(spatial_data(mixed_data, region))
        @test result isa AbstractMatrix
        @test size(result) == (3, 3)
        @test isdiag(result)

        # First should be point process (positive)
        @test result[1, 1] > 0
        # Second should be grid (zero)
        @test result[2, 2] == 0.0
        # Third should be marked points (positive)
        @test result[3, 3] > 0
    end

    @testset "One dimensional test" begin
        # Test different dimensions
        points_1d, region_1d = make_points_example(rng, n_processes = 1,
            region_min = (-2.0,),
            region_max = (3.0,),
            point_number = 20)

        result_1d = covariance_zero_atom(spatial_data(points_1d[1], region_1d))
        @test isfinite(result_1d)
        @test result_1d > 0
    end
end
