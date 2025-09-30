using Test, SpatialMultitaper, StableRNGs
include("../test_utilities/TestData.jl")
using .TestData

import SpatialMultitaper: _validate_nk, _validate_dk, _validate_kmax, _validate_dims,
                          _validate_wavenumber_params, default_dk

@testset "Parameter Input Tests" begin
    @testset "_validate_dims tests" begin
        @testset "Tuple input" begin
            @test _validate_dims((1, 2, 3), 3) == (1, 2, 3)
            @test _validate_dims((0.1, 0.2), 2) == (0.1, 0.2)

            # Wrong dimension
            @test_throws ArgumentError _validate_dims((1, 2), 3)
            @test_throws ArgumentError _validate_dims((1, 2, 3, 4), 2)
        end

        @testset "Number input" begin
            @test _validate_dims(5, 3) == (5, 5, 5)
            @test _validate_dims(0.1, 2) == (0.1, 0.1)
            @test _validate_dims(42, 1) == (42,)
        end
    end

    @testset "_validate_nk tests" begin
        @testset "Valid inputs" begin
            @test _validate_nk((32, 64), 2) == (32, 64)
            @test _validate_nk(128, 3) == (128, 128, 128)
            @test _validate_nk((16, 32, 64), 3) == (16, 32, 64)
        end

        @testset "Invalid inputs" begin
            # Non-integer
            @test_throws ArgumentError _validate_nk((32.5, 64), 2)
            @test_throws ArgumentError _validate_nk(32.0, 2)

            # Non-positive
            @test_throws ArgumentError _validate_nk((0, 64), 2)
            @test_throws ArgumentError _validate_nk((-5, 32), 2)
            @test_throws ArgumentError _validate_nk(-10, 2)

            # Wrong dimension
            @test_throws ArgumentError _validate_nk((32, 64), 3)
        end
    end

    @testset "_validate_dk tests" begin
        @testset "Valid inputs" begin
            @test _validate_dk((0.1, 0.2), 2) == (0.1, 0.2)
            @test _validate_dk(0.05, 3) == (0.05, 0.05, 0.05)
            @test _validate_dk((0.01, 0.02, 0.03), 3) == (0.01, 0.02, 0.03)
        end

        @testset "Invalid inputs" begin
            # Non-positive
            @test_throws ArgumentError _validate_dk((0.0, 0.1), 2)
            @test_throws ArgumentError _validate_dk((-0.1, 0.2), 2)
            @test_throws ArgumentError _validate_dk(-0.05, 2)

            # Wrong dimension
            @test_throws ArgumentError _validate_dk((0.1, 0.2), 3)
        end
    end

    @testset "_validate_kmax tests" begin
        @testset "Valid inputs" begin
            @test _validate_kmax((1.0, 2.0), 2) == (1.0, 2.0)
            @test _validate_kmax(0.5, 3) == (0.5, 0.5, 0.5)
            @test _validate_kmax((0.1, 0.2, 0.3), 3) == (0.1, 0.2, 0.3)
        end

        @testset "Invalid inputs" begin
            # Non-positive
            @test_throws ArgumentError _validate_kmax((0.0, 1.0), 2)
            @test_throws ArgumentError _validate_kmax((-0.5, 1.0), 2)
            @test_throws ArgumentError _validate_kmax(-1.0, 2)

            # Wrong dimension
            @test_throws ArgumentError _validate_kmax((1.0, 2.0), 3)
        end
    end

    @testset "_validate_wavenumber_params tests" begin
        @testset "nk and kmax provided" begin
            nk, kmax = _validate_wavenumber_params((32, 64), (1.0, 2.0), nothing, 2)
            @test nk == (32, 64)
            @test kmax == (1.0, 2.0)

            # Single values
            nk, kmax = _validate_wavenumber_params(128, 0.5, nothing, 2)
            @test nk == (128, 128)
            @test kmax == (0.5, 0.5)
        end

        @testset "nk and dk provided" begin
            nk, kmax = _validate_wavenumber_params((32, 64), nothing, (0.1, 0.2), 2)
            @test nk == (32, 64)
            @test kmax == (32 * 0.1 / 2, 64 * 0.2 / 2)
            @test kmax == (1.6, 6.4)

            # Single values
            nk, kmax = _validate_wavenumber_params(100, nothing, 0.02, 2)
            @test nk == (100, 100)
            @test kmax == (1.0, 1.0)
        end

        @testset "kmax and dk provided" begin
            nk, kmax = _validate_wavenumber_params(nothing, (1.0, 2.0), (0.1, 0.2), 2)
            @test kmax == (1.0, 2.0)
            # nk = ceil(2 * kmax / dk)
            @test nk == (ceil(2 * 1.0 / 0.1), ceil(2 * 2.0 / 0.2))
            @test nk == (20, 20)
            @test eltype(nk) == Int

            # Single values
            nk, kmax = _validate_wavenumber_params(nothing, 0.5, 0.05, 2)
            @test kmax == (0.5, 0.5)
            @test nk == (ceil(2 * 0.5 / 0.05), ceil(2 * 0.5 / 0.05))
            @test nk == (20, 20)
            @test eltype(nk) == Int
        end

        @testset "Error cases - insufficient parameters" begin
            # No parameters
            @test_throws ErrorException _validate_wavenumber_params(
                nothing, nothing, nothing, 2)

            # Only one parameter
            @test_throws ErrorException _validate_wavenumber_params(32, nothing, nothing, 2)
            @test_throws ErrorException _validate_wavenumber_params(
                nothing, 1.0, nothing, 2)
            @test_throws ErrorException _validate_wavenumber_params(
                nothing, nothing, 0.1, 2)
        end

        @testset "Consistency checks" begin
            # When all three are mathematically consistent
            nk, kmax = _validate_wavenumber_params((40, 80), (1.0, 2.0), nothing, 2)
            dk_calc = 2 .* kmax ./ nk
            @test all(dk_calc .≈ (0.05, 0.05))

            # Verify relationship: kmax = nk * dk / 2
            nk_test = (50, 100)
            dk_test = (0.04, 0.02)
            nk_out, kmax_out = _validate_wavenumber_params(nk_test, nothing, dk_test, 2)
            @test nk_out == nk_test
            @test all(kmax_out .≈ nk_test .* dk_test ./ 2)
        end
    end

    @testset "default_dk tests" begin
        @testset "From region" begin
            # Create a simple rectangular region
            rng = StableRNG(123)
            points_data = make_points_example(rng, n_processes = 1, dim = 2,
                return_type = :single, point_number = 10)
            region = getregion(points_data)

            dk_default = default_dk(region, nothing, nothing)
            @test dk_default isa NTuple{2}
            @test all(dk_default .> 0)
        end

        @testset "From SpatialData" begin
            rng = StableRNG(123)
            spatial_data = make_points_example(rng, n_processes = 1, dim = 2,
                return_type = :single, point_number = 20)

            dk_default = default_dk(spatial_data, nothing, nothing)
            @test dk_default isa NTuple{2}
            @test all(dk_default .> 0)
        end

        @testset "Different dimensions" begin
            rng = StableRNG(123)

            # 1D case
            data_1d = make_points_example(rng, n_processes = 1, dim = 1,
                return_type = :single, point_number = 15)
            dk_1d = default_dk(data_1d, nothing, nothing)
            @test dk_1d isa NTuple{1}
            @test dk_1d[1] > 0

            # 3D case
            data_3d = make_points_example(rng, n_processes = 1, dim = 3,
                return_type = :single, point_number = 15)
            dk_3d = default_dk(data_3d, nothing, nothing)
            @test dk_3d isa NTuple{3}
            @test all(dk_3d .> 0)

            @test isnothing(default_dk(data_3d, 1, 1.0))
        end
    end

    @testset "Edge cases and type stability" begin
        @testset "Large numbers" begin
            @test _validate_nk(1000000, 1) == (1000000,)
            @test _validate_kmax(1e6, 1) == (1e6,)
            @test _validate_dk(1e-10, 1) == (1e-10,)
        end

        @testset "Type preservation" begin
            # Integer types
            @test typeof(_validate_nk(Int32(64), 1)[1]) == Int32
            @test typeof(_validate_nk(Int64(64), 1)[1]) == Int64

            # Float types
            @test typeof(_validate_dk(Float32(0.1), 1)[1]) == Float32
            @test typeof(_validate_kmax(Float64(1.0), 1)[1]) == Float64
        end

        @testset "Mixed dimension validation" begin
            # Test that all validation functions handle dimensions correctly
            for dim in 1:5
                @test length(_validate_nk(64, dim)) == dim
                @test length(_validate_dk(0.1, dim)) == dim
                @test length(_validate_kmax(1.0, dim)) == dim
            end
        end
    end
end
