using Test, SpatialMultitaper, StaticArrays, LinearAlgebra, StableRNGs

import SpatialMultitaper: covariance_zero_atom, unitless_measure

# Import test utilities
include("test_utilities/TestData.jl")
using .TestData

@testset "covariance_zero_atom tests" begin
    rng = StableRNG(123)

    @testset "PointSet tests" begin
        points_data = make_points_example(
            rng, n_processes = 1, return_type = :single, point_number = 100)

        # Test basic computation
        result = covariance_zero_atom(points_data)
        expected = length(observations(points_data)) /
                   unitless_measure(getregion(points_data))
        @test result ≈ expected
        @test result > 0  # Should be positive for non-empty point set

        # Test with different point densities
        sparse_data = make_points_example(
            rng, n_processes = 1, return_type = :single, point_number = 10)
        sparse_result = covariance_zero_atom(sparse_data)
        @test sparse_result < result  # Lower density should give smaller atom
    end

    @testset "Marked point tests (GeoTable)" begin
        marked_data = make_marked_example(
            rng, n_processes = 1, return_type = :single, point_number = 50)

        result = covariance_zero_atom(marked_data)

        # Should compute sum of squared marks divided by region measure
        marks = values(observations(marked_data))[1]  # Get the marks
        expected = sum(abs2, marks) / unitless_measure(getregion(marked_data))
        @test result ≈ expected
        @test result > 0  # Should be positive for non-zero marks
    end

    @testset "CartesianGrid tests" begin
        grids_data = make_grids_example(
            rng, n_processes = 1, return_type = :single, grid_dims = (10, 8))

        # Grid data should always return zero atom
        result = covariance_zero_atom(grids_data)
        @test result == 0.0

        # Test with different grid sizes
        large_grids = make_grids_example(
            rng, n_processes = 1, return_type = :single, grid_dims = (50, 30))
        large_result = covariance_zero_atom(large_grids)
        @test large_result == 0.0
    end

    @testset "Multiple processes - Tuple format" begin
        # Test with multiple point processes
        points_data = make_points_example(
            rng, n_processes = 3, return_type = :tuple, point_number = 60)

        result = covariance_zero_atom(points_data)
        @test result isa AbstractMatrix
        @test size(result) == (3, 3)
        @test isdiag(result)  # Should be diagonal matrix

        # Diagonal elements should be individual atom estimates
        for i in 1:3
            individual_data = make_points_example(
                rng, n_processes = 1, return_type = :single, point_number = 60)
            individual_result = covariance_zero_atom(individual_data)
            @test result[i, i] > 0
        end
    end

    @testset "Multiple processes - Vector format" begin
        points_data = make_points_example(rng, n_processes = 2, return_type = :vector)

        result = covariance_zero_atom(points_data)
        @test result isa AbstractMatrix
        @test size(result) == (2, 2)
        @test isdiag(result)

        # Check that diagonal elements are positive
        @test result[1, 1] > 0
        @test result[2, 2] > 0
    end

    @testset "Single element tuple" begin
        points_data = make_points_example(rng, n_processes = 1, return_type = :tuple)

        result = covariance_zero_atom(points_data)
        expected_single = make_points_example(rng, n_processes = 1, return_type = :single)
        expected = covariance_zero_atom(expected_single)
        @test result isa SMatrix
        @test result[1, 1] > 0
    end

    @testset "Mixed data types" begin
        mixed_data = make_mixed_example(rng, n_processes = (1, 1, 1), return_type = :tuple)

        result = covariance_zero_atom(mixed_data)
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
        points_1d = make_points_example(rng, n_processes = 1, return_type = :single,
            dim = 1, region_min = (-2.0,), region_max = (3.0,), point_number = 20)

        result_1d = covariance_zero_atom(points_1d)
        @test isfinite(result_1d)
        @test result_1d > 0
    end
end
