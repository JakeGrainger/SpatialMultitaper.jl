using Test, SpatialMultitaper, StableRNGs, StaticArrays, LinearAlgebra, BenchmarkTools
include("../../test_utilities/TestData.jl")
using .TestData

import SpatialMultitaper: getestimates, getevaluationpoints, KFunction, k_function!,
                          getprocessinformation, process_trait, _c_to_k_transform!
rng = StableRNG(123)

# loop over 1d, 2d, 3d
@testset "Dimension $dim tests" for dim in 1:3

    # - k function from raw data is same as from Spatial data
    @testset "Raw data vs SpatialData consistency" begin
        # Test point processes
        points_raw, region = make_points_example(rng, n_processes = 1, dim = dim,
            return_type = :raw, point_number = 50)
        points_spatial = spatial_data(points_raw, region)

        # Define test parameters
        radii = dim == 1 ? [0.1, 0.5, 1.0] :
                dim == 2 ? [0.1, 0.3, 0.8] : [0.1, 0.2, 0.5]
        nk = dim == 1 ? (32,) : dim == 2 ? (16, 16) : (8, 8, 8)
        kmax = dim == 1 ? (2.0,) : dim == 2 ? (1.5, 1.5) : (1.0, 1.0, 1.0)
        bw = ntuple(_ -> 3, dim)
        tapers = sin_taper_family(bw, region)

        raw_result = k_function(
            points_raw, region, radii = radii, nk = nk, kmax = kmax, tapers = tapers)
        spatial_result = k_function(
            points_spatial, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

        @test getevaluationpoints(raw_result) ≈ getevaluationpoints(spatial_result)
        @test getestimates(raw_result) ≈ getestimates(spatial_result)
    end

    # - loop over single, tuple and vector
    @testset "Return type: $return_type" for return_type in [:single, :tuple, :vector]
        n_processes = return_type == :single ? 1 : 3

        # - - k function from SpatialData
        @testset "K-function from SpatialData" begin
            # Point processes
            points_data = make_points_example(
                rng, n_processes = n_processes, dim = dim,
                return_type = return_type, point_number = 40)

            radii = dim == 1 ? [0.2, 0.6] : dim == 2 ? [0.2, 0.5] : [0.1, 0.3]
            nk = dim == 1 ? (24,) : dim == 2 ? (12, 12) : (6, 6, 6)
            kmax = dim == 1 ? (1.5,) : dim == 2 ? (1.2, 1.2) : (0.8, 0.8, 0.8)
            bw = ntuple(_ -> 3, dim)
            region = getregion(points_data)
            tapers = sin_taper_family(bw, region)

            result = k_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test result isa KFunction
            @test getevaluationpoints(result) ≈ radii
            if return_type == :vector
                @test size(getestimates(result)) ==
                      (n_processes, n_processes, length(radii))
            else
                @test length(getestimates(result)) == length(radii)
            end
            @test all(x -> all(isfinite.(x)), getestimates(result))
        end

        # - - k function from spectra
        @testset "K-function from spectra" begin
            points_data = make_points_example(
                rng, n_processes = n_processes, dim = dim,
                return_type = return_type, point_number = 35)

            radii = dim == 1 ? [0.15, 0.4] : [0.15, 0.35]
            nk = dim == 1 ? (20,) : dim == 2 ? (10, 10) : (6, 6, 6)
            kmax = dim == 1 ? (1.2,) : dim == 2 ? (1.0, 1.0) : (0.6, 0.6, 0.6)
            bw = ntuple(_ -> 3, dim)
            region = getregion(points_data)
            tapers = sin_taper_family(bw, region)

            # First compute spectra
            spectrum = spectra(points_data, nk = nk, kmax = kmax, tapers = tapers)

            # Then compute k_function from spectra
            k_from_spectra = k_function(spectrum, radii = radii)

            # Compare with direct computation
            k_direct = k_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test getevaluationpoints(k_from_spectra) ≈ getevaluationpoints(k_direct)
            @test getestimates(k_from_spectra) ≈ getestimates(k_direct)
        end

        # - - k function from c function
        @testset "K-function from C-function" begin
            points_data = make_points_example(
                rng, n_processes = n_processes, dim = dim,
                return_type = return_type, point_number = 35)

            radii = dim == 1 ? [0.15, 0.4] : [0.15, 0.35]
            nk = dim == 1 ? (20,) : dim == 2 ? (10, 10) : (6, 6, 6)
            kmax = dim == 1 ? (1.2,) : dim == 2 ? (1.0, 1.0) : (0.6, 0.6, 0.6)
            bw = ntuple(_ -> 3, dim)
            region = getregion(points_data)
            tapers = sin_taper_family(bw, region)

            # First compute C function
            c_func = c_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            # Then compute K function from C function
            k_from_c = k_function(c_func)

            # Compare with direct computation
            k_direct = k_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test getevaluationpoints(k_from_c) ≈ getevaluationpoints(k_direct)
            @test getestimates(k_from_c) ≈ getestimates(k_direct)
        end

        # - - partial k function from SpatialData
        @testset "Partial k-function from SpatialData" begin
            if n_processes > 1  # Partial only makes sense for multiple processes
                points_data = make_points_example(
                    rng, n_processes = n_processes, dim = dim,
                    return_type = return_type, point_number = 30)

                radii = dim == 1 ? [0.1, 0.3] : [0.1, 0.25]
                nk = dim == 1 ? (16,) : dim == 2 ? (8, 8) : (4, 4, 4)
                kmax = dim == 1 ? (1.0,) : dim == 2 ? (0.8, 0.8) : (0.5, 0.5, 0.5)
                bw = ntuple(_ -> 3, dim)
                region = getregion(points_data)
                tapers = sin_taper_family(bw, region)

                partial_result = partial_k_function(points_data, radii = radii,
                    nk = nk, kmax = kmax, tapers = tapers)

                @test partial_result isa KFunction
                @test getevaluationpoints(partial_result) ≈ radii
                if return_type == :vector
                    @test size(getestimates(partial_result)) ==
                          (n_processes, n_processes, length(radii))
                else
                    @test length(getestimates(partial_result)) == length(radii)
                end
            end
        end

        # - - partial k function from partial spectra
        @testset "Partial k-function from partial spectra" begin
            if n_processes > 1
                points_data = make_points_example(
                    rng, n_processes = n_processes, dim = dim,
                    return_type = return_type, point_number = 25)

                radii = dim == 1 ? [0.2] : [0.2]
                nk = dim == 1 ? (12,) : dim == 2 ? (8, 8) : (4, 4, 4)
                kmax = dim == 1 ? (0.8,) : dim == 2 ? (0.6, 0.6) : (0.4, 0.4, 0.4)
                bw = ntuple(_ -> 3, dim)
                region = getregion(points_data)
                tapers = sin_taper_family(bw, region)

                # Compute partial spectra first
                partial_spec = partial_spectra(
                    points_data, nk = nk, kmax = kmax, tapers = tapers)

                # Then k_function from partial spectra
                k_from_partial = k_function(partial_spec, radii = radii)

                @test k_from_partial isa KFunction
                @test getevaluationpoints(k_from_partial) ≈ radii
            end
        end

        # - - partial k function from partial c function
        @testset "Partial k-function from partial C-function" begin
            if n_processes > 1
                points_data = make_points_example(
                    rng, n_processes = n_processes, dim = dim,
                    return_type = return_type, point_number = 25)

                radii = dim == 1 ? [0.2] : [0.2]
                nk = dim == 1 ? (12,) : dim == 2 ? (8, 8) : (4, 4, 4)
                kmax = dim == 1 ? (0.8,) : dim == 2 ? (0.6, 0.6) : (0.4, 0.4, 0.4)
                bw = ntuple(_ -> 3, dim)
                region = getregion(points_data)
                tapers = sin_taper_family(bw, region)

                # Compute partial C function first
                partial_c = partial_c_function(
                    points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

                # Then K function from partial C function
                k_from_partial_c = k_function(partial_c)

                @test k_from_partial_c isa KFunction
                @test getevaluationpoints(k_from_partial_c) ≈ radii

                # Compare with direct partial K function computation
                k_partial_direct = partial_k_function(
                    points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

                @test getevaluationpoints(k_from_partial_c) ≈
                      getevaluationpoints(k_partial_direct)
                @test getestimates(k_from_partial_c) ≈ getestimates(k_partial_direct)
            end
        end

        # - - check stored values match the types correctly, so Number, SMatrix, Array of Number with dim(array) = D+2
        @testset "Type correctness" begin
            points_data = make_points_example(
                rng, n_processes = n_processes, dim = dim,
                return_type = return_type, point_number = 25)

            radii = [0.1, 0.2, 0.3]
            nk = dim == 1 ? (8,) : dim == 2 ? (6, 6) : (4, 4, 4)
            kmax = dim == 1 ? (0.5,) : dim == 2 ? (0.4, 0.4) : (0.3, 0.3, 0.3)
            bw = ntuple(_ -> 3, dim)
            region = getregion(points_data)
            tapers = sin_taper_family(bw, region)

            result = k_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test result isa KFunction
            @test getevaluationpoints(result) isa AbstractVector
            @test length(getevaluationpoints(result)) == length(radii)
            if return_type == :vector
                @test getestimates(result) isa AbstractArray
                @test size(getestimates(result)) ==
                      (n_processes, n_processes, length(radii))
            elseif return_type == :tuple
                @test getestimates(result) isa AbstractVector
                @test length(getestimates(result)) == length(radii)
                @test eltype(getestimates(result)) <: SMatrix{n_processes, n_processes}
            else # return_type == :single
                @test getestimates(result) isa AbstractVector
                @test length(getestimates(result)) == length(radii)
            end

            # Check element types
            @test eltype(getevaluationpoints(result)) <: Number
            if return_type == :tuple
                @test eltype(getestimates(result)) <: SMatrix
            else
                @test eltype(getestimates(result)) <: Number
            end
        end

        # - - check indexing works
        @testset "Indexing functionality" begin
            points_data = make_points_example(
                rng, n_processes = n_processes, dim = dim,
                return_type = return_type, point_number = 30)

            radii = [0.1, 0.25, 0.5]
            nk = dim == 1 ? (12,) : dim == 2 ? (8, 8) : (4, 4, 4)
            kmax = dim == 1 ? (0.8,) : dim == 2 ? (0.6, 0.6) : (0.4, 0.4, 0.4)
            bw = ntuple(_ -> 3, dim)
            region = getregion(points_data)
            tapers = sin_taper_family(bw, region)

            result = k_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            # Test indexing into results
            radii_result = getevaluationpoints(result)
            values_result = getestimates(result)

            @test radii_result[1] ≈ radii[1]
            @test radii_result[end] ≈ radii[end]
            @test length(radii_result) == length(radii)

            if return_type == :tuple
                @test values_result[1] isa SMatrix{n_processes, n_processes}
                @test values_result[end] isa SMatrix{n_processes, n_processes}
            elseif return_type == :vector
                @test size(values_result) == (n_processes, n_processes, length(radii))
            else # return_type == :single
                @test values_result[1] isa Number
                @test values_result[end] isa Number
            end

            if return_type == :vector
                @test size(values_result) == (n_processes, n_processes, length(radii))
            else
                @test length(values_result) == length(radii)
            end

            # Test that values are finite
            @test all(x -> all(isfinite.(x)), values_result)
        end

        # - - K function specific tests (relationship to C function)
        @testset "K-function specific properties" begin
            points_data = make_points_example(
                rng, n_processes = n_processes, dim = dim,
                return_type = return_type, point_number = 30)

            radii = dim == 1 ? [0.1, 0.3] : dim == 2 ? [0.1, 0.3] : [0.1, 0.2]
            nk = dim == 1 ? (16,) : dim == 2 ? (10, 10) : (6, 6, 6)
            kmax = dim == 1 ? (1.0,) : dim == 2 ? (0.8, 0.8) : (0.6, 0.6, 0.6)
            bw = ntuple(_ -> 3, dim)
            region = getregion(points_data)
            tapers = sin_taper_family(bw, region)

            # Compute both K and C functions
            k_result = k_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)
            c_result = c_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test getevaluationpoints(k_result) ≈ getevaluationpoints(c_result)

            # Test that K function computed from C function gives same result
            k_from_c = k_function(c_result)
            @test getevaluationpoints(k_from_c) ≈ getevaluationpoints(k_result)
            @test getestimates(k_from_c) ≈ getestimates(k_result)

            # K function should be non-negative for positive radii and single process
            if return_type == :single
                k_values = getestimates(k_result)
                @test all(k_values .>= 0) ||
                      length(findall(x -> x < 0, k_values)) / length(k_values) < 0.1  # Allow small numerical errors
            end
        end
    end
end

@testset "Edge cases and robustness" begin
    # Test very small radii
    points_data = make_points_example(rng, n_processes = 1, dim = 2,
        return_type = :single, point_number = 20)
    small_radii = [0.01, 0.05]
    bw = (3, 3)
    region = getregion(points_data)
    tapers = sin_taper_family(bw, region)
    result = k_function(
        points_data, radii = small_radii, nk = (8, 8), kmax = (0.5, 0.5), tapers = tapers)
    @test all(x -> all(isfinite.(x)), getestimates(result))

    # Test large radii
    large_radii = [2.0, 100.0]
    @test_throws ArgumentError k_function(
        points_data, radii = large_radii, nk = (8, 8), kmax = (0.5, 0.5), tapers = tapers)

    # Test K function at zero (should be zero or near zero)
    zero_radius = [0.0]
    result_zero = k_function(
        points_data, radii = zero_radius, nk = (8, 8), kmax = (0.5, 0.5), tapers = tapers)
    @test abs(getestimates(result_zero)[1]) < 1e-6  # K(0) should be approximately 0
end

@testset "in place tests $return_type" for return_type in [:single, :tuple, :vector]
    points_data = make_points_example(rng, n_processes = 1, dim = 2,
        return_type = :single, point_number = 20)
    tapers = sin_taper_family((3, 3), getregion(points_data))
    radii = [0.1, 0.3, 0.5]
    spectrum = spectra(points_data, nk = (16, 16), kmax = (0.5, 0.5), tapers = tapers)
    c_fun = c_function(spectrum, radii = radii)

    mean_prod = getprocessinformation(c_fun).mean_product
    alloc = @ballocated _c_to_k_transform!(
        $radii, getestimates($c_fun), process_trait($spectrum), $mean_prod, Val{2}()) samples=1
    @test alloc == 0

    # from spatial data
    alloc = @ballocated k_function!($spectrum, radii = $radii) samples=1
    @test_broken alloc == 0

    # from c function
    alloc = @ballocated k_function!($c_fun) samples=1
    @test alloc == 0
end
