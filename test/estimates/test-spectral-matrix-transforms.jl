using SpatialMultitaper, Test, StableRNGs, StaticArrays, LinearAlgebra
include("../test_utilities/TestData.jl")
using .TestData

import SpatialMultitaper: apply_transform, is_partial, Spectra, Coherence

@testset "apply_transform function" begin
    @testset "SMatrix input" begin
        # Create test data with SMatrix elements
        data_xy = (1:4, 1:4)
        data = [SMatrix{2, 2}(rand(ComplexF64, 2, 2)) for _ in 1:4, _ in 1:4]

        @testset "Function with no additional args" begin
            result = apply_transform(real, data_xy, data)
            @test size(result) == size(data)
            @test eltype(result) == SMatrix{2, 2, Float64, 4}
            @test all(x -> all(imag.(x) .≈ 0), result)
        end

        @testset "Function with additional args" begin
            # Mock function that uses additional arguments
            mock_transform(x, factor) = factor .* x
            result = apply_transform(mock_transform, data_xy, data, 2.0)
            @test size(result) == size(data)
            @test all(x ≈ 2.0 * y for (x, y) in zip(result, data))
        end
    end

    @testset "Regular array input" begin
        data_xy = (1:6, 1:6)
        data = rand(ComplexF64, 3, 3, 6, 6)  # P x Q x freq1 x freq2

        @testset "Transform along matrix dimensions" begin
            # Test with function that operates on matrices
            det_transform(x) = det(x)
            result = apply_transform(det_transform, data_xy, data)

            @test size(result) == (1, 1, 6, 6)  # Should not collapse first two dimensions
            @test eltype(result) <: ComplexF64
        end

        @testset "Transform with args" begin
            matrix_power(x, p) = x^p
            result = apply_transform(matrix_power, data_xy, data, 2)

            @test size(result) == (3, 3, 6, 6)
        end
    end

    @testset "Edge cases" begin
        # Single element arrays
        single_element = [SMatrix{1, 1}(rand(ComplexF64, 1, 1))]
        result = apply_transform(x -> abs.(x), (1:1,), single_element)
        @test length(result) == 1
        @test eltype(result) == SMatrix{1, 1, Float64, 1}

        # Empty arrays
        empty_data = SMatrix{2, 2, ComplexF64, 4}[]
        result = apply_transform(real, (Int[],), empty_data)
        @test isempty(result)
    end
end

@testset "Integration with Estimate Types" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)

    @testset "apply_transform through coherence calculations" begin
        region = getregion(data)
        spec = spectra(data, nfreq = (4, 4), fmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        # Test that coherence uses apply_transform internally
        coh = coherence(spec)
        @test coh isa Coherence
        @test size(coh) == size(spec)
    end

    @testset "apply_transform through partial_spectra" begin
        region = getregion(data)
        spec = spectra(data, nfreq = (4, 4), fmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        # Test that partial_spectra uses apply_transform internally
        partial_spec = partial_spectra(spec)
        @test partial_spec isa Spectra
        @test is_partial(partial_spec) == true
    end
end

@testset "Type Preservation" begin
    @testset "Float32 preservation" begin
        data_xy = (1:3, 1:3)
        data = [SMatrix{2, 2}(rand(Float32, 2, 2)) for _ in data_xy[1], _ in data_xy[2]]
        result = apply_transform(x -> 2 .* x, data_xy, data)
        @test eltype(eltype(result)) == Float32
    end

    @testset "Complex type handling" begin
        data_xy = (1:3, 1:3)
        data = [SMatrix{2, 2}(rand(ComplexF32, 2, 2)) for _ in data_xy[1], _ in data_xy[2]]
        result_real = apply_transform(x -> real.(x), data_xy, data)
        @test eltype(eltype(result_real)) == Float32

        result_abs = apply_transform(x -> abs.(x), data_xy, data)
        @test eltype(eltype(result_abs)) == Float32
    end
end
