using Test, SpatialMultitaper, StableRNGs, StaticArrays, LinearAlgebra, BenchmarkTools
include("../../test_utilities/TestData.jl")
using .TestData

import SpatialMultitaper: getestimate, getevaluationpoints, CenteredLFunction,
                          centered_l_function!

#
rng = StableRNG(123)

# loop over 1d, 2d, 3d
@testset "Dimension $dim tests" for dim in 1:3

    # - centered l function from raw data is same as from Spatial data
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

        raw_result = centered_l_function(
            points_raw, region, radii = radii, nk = nk, kmax = kmax, tapers = tapers)
        spatial_result = centered_l_function(
            points_spatial, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

        @test getevaluationpoints(raw_result) ≈ getevaluationpoints(spatial_result)
        @test getestimate(raw_result) ≈ getestimate(spatial_result)
    end

    # - loop over single, tuple and vector
    @testset "Return type: $return_type" for return_type in [:single, :tuple, :vector]
        n_processes = return_type == :single ? 1 : 3

        # - - centered l function from SpatialData
        @testset "Centered L-function from SpatialData" begin
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

            result = centered_l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test result isa CenteredLFunction
            @test getevaluationpoints(result) ≈ radii
            if return_type == :vector
                @test size(getestimate(result)) ==
                      (n_processes, n_processes, length(radii))
            else
                @test length(getestimate(result)) == length(radii)
            end
            @test all(x -> all(isfinite.(x)), getestimate(result))
        end

        # - - centered l function from spectra
        @testset "Centered L-function from spectra" begin
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

            # Then compute centered_l_function from spectra
            centered_l_from_spectra = centered_l_function(spectrum, radii = radii)

            # Compare with direct computation
            centered_l_direct = centered_l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test getevaluationpoints(centered_l_from_spectra) ≈
                  getevaluationpoints(centered_l_direct)
            @test getestimate(centered_l_from_spectra) ≈ getestimate(centered_l_direct)
        end

        # - - centered l function from c function
        @testset "Centered L-function from C-function" begin
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

            # Then compute centered L function from C function
            centered_l_from_c = centered_l_function(c_func)

            # Compare with direct computation
            centered_l_direct = centered_l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test getevaluationpoints(centered_l_from_c) ≈
                  getevaluationpoints(centered_l_direct)
            @test getestimate(centered_l_from_c) ≈ getestimate(centered_l_direct)
        end

        # - - centered l function from k function
        @testset "Centered L-function from K-function" begin
            points_data = make_points_example(
                rng, n_processes = n_processes, dim = dim,
                return_type = return_type, point_number = 35)

            radii = dim == 1 ? [0.15, 0.4] : [0.15, 0.35]
            nk = dim == 1 ? (20,) : dim == 2 ? (10, 10) : (6, 6, 6)
            kmax = dim == 1 ? (1.2,) : dim == 2 ? (1.0, 1.0) : (0.6, 0.6, 0.6)
            bw = ntuple(_ -> 3, dim)
            region = getregion(points_data)
            tapers = sin_taper_family(bw, region)

            # First compute K function
            k_func = k_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            # Then compute centered L function from K function
            centered_l_from_k = centered_l_function(k_func)

            # Compare with direct computation
            centered_l_direct = centered_l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test getevaluationpoints(centered_l_from_k) ≈
                  getevaluationpoints(centered_l_direct)
            @test getestimate(centered_l_from_k) ≈ getestimate(centered_l_direct)
        end

        # - - centered l function from l function
        @testset "Centered L-function from L-function" begin
            points_data = make_points_example(
                rng, n_processes = n_processes, dim = dim,
                return_type = return_type, point_number = 35)

            radii = dim == 1 ? [0.15, 0.4] : [0.15, 0.35]
            nk = dim == 1 ? (20,) : dim == 2 ? (10, 10) : (6, 6, 6)
            kmax = dim == 1 ? (1.2,) : dim == 2 ? (1.0, 1.0) : (0.6, 0.6, 0.6)
            bw = ntuple(_ -> 3, dim)
            region = getregion(points_data)
            tapers = sin_taper_family(bw, region)

            # First compute L function
            l_func = l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            # Then compute centered L function from L function
            centered_l_from_l = centered_l_function(l_func)

            # Compare with direct computation
            centered_l_direct = centered_l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test getevaluationpoints(centered_l_from_l) ≈
                  getevaluationpoints(centered_l_direct)
            @test getestimate(centered_l_from_l) ≈ getestimate(centered_l_direct)
        end

        # - - partial centered l function from SpatialData
        @testset "Partial centered L-function from SpatialData" begin
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

                partial_result = partial_centered_l_function(
                    points_data, radii = radii,
                    nk = nk, kmax = kmax, tapers = tapers)

                @test partial_result isa CenteredLFunction
                @test getevaluationpoints(partial_result) ≈ radii
                if return_type == :vector
                    @test size(getestimate(partial_result)) ==
                          (n_processes, n_processes, length(radii))
                else
                    @test length(getestimate(partial_result)) == length(radii)
                end
            end
        end

        # - - partial centered l function from partial functions
        @testset "Partial centered L-function from partial functions" begin
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

                # Compute partial functions first
                partial_c = partial_c_function(
                    points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)
                partial_k = partial_k_function(
                    points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)
                partial_l = partial_l_function(
                    points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

                # Then centered L functions from each
                centered_l_from_partial_c = centered_l_function(partial_c)
                centered_l_from_partial_k = centered_l_function(partial_k)
                centered_l_from_partial_l = centered_l_function(partial_l)

                @test centered_l_from_partial_c isa CenteredLFunction
                @test centered_l_from_partial_k isa CenteredLFunction
                @test centered_l_from_partial_l isa CenteredLFunction

                # All should give same result
                @test getevaluationpoints(centered_l_from_partial_c) ≈
                      getevaluationpoints(centered_l_from_partial_k)
                @test getevaluationpoints(centered_l_from_partial_c) ≈
                      getevaluationpoints(centered_l_from_partial_l)
                @test getestimate(centered_l_from_partial_c) ≈
                      getestimate(centered_l_from_partial_k)
                @test getestimate(centered_l_from_partial_c) ≈
                      getestimate(centered_l_from_partial_l)
            end
        end

        # - - check stored values match the types correctly
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

            result = centered_l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test result isa CenteredLFunction
            @test getevaluationpoints(result) isa AbstractVector
            @test length(getevaluationpoints(result)) == length(radii)
            if return_type == :vector
                @test getestimate(result) isa AbstractArray
                @test size(getestimate(result)) ==
                      (n_processes, n_processes, length(radii))
            elseif return_type == :tuple
                @test getestimate(result) isa AbstractVector
                @test length(getestimate(result)) == length(radii)
                @test eltype(getestimate(result)) <: SMatrix{n_processes, n_processes}
            else # return_type == :single
                @test getestimate(result) isa AbstractVector
                @test length(getestimate(result)) == length(radii)
            end

            # Check element types
            @test eltype(getevaluationpoints(result)) <: Number
            if return_type == :tuple
                @test eltype(getestimate(result)) <: SMatrix
            else
                @test eltype(getestimate(result)) <: Number
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

            result = centered_l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            # Test indexing into results
            radii_result = getevaluationpoints(result)
            values_result = getestimate(result)

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

        # - - Centered L function specific tests (relationship to L function)
        @testset "Centered L-function specific properties" begin
            points_data = make_points_example(
                rng, n_processes = n_processes, dim = dim,
                return_type = return_type, point_number = 30)

            radii = dim == 1 ? [0.1, 0.3] : dim == 2 ? [0.1, 0.3] : [0.1, 0.2]
            nk = dim == 1 ? (16,) : dim == 2 ? (10, 10) : (6, 6, 6)
            kmax = dim == 1 ? (1.0,) : dim == 2 ? (0.8, 0.8) : (0.6, 0.6, 0.6)
            bw = ntuple(_ -> 3, dim)
            region = getregion(points_data)
            tapers = sin_taper_family(bw, region)

            # Compute centered L and L functions
            centered_l_result = centered_l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)
            l_result = l_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test getevaluationpoints(centered_l_result) ≈ getevaluationpoints(l_result)

            # Test that centered L function computed from L function gives same result
            centered_l_from_l = centered_l_function(l_result)
            @test getevaluationpoints(centered_l_from_l) ≈
                  getevaluationpoints(centered_l_result)
            @test getestimate(centered_l_from_l) ≈ getestimate(centered_l_result)

            # Verify the centering relationship: centered_L = L - r
            if return_type == :single
                l_values = getestimate(l_result)
                centered_l_values = getestimate(centered_l_result)
                expected_centered = l_values .- radii
                @test centered_l_values ≈ expected_centered
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
    result = centered_l_function(
        points_data, radii = small_radii, nk = (8, 8), kmax = (0.5, 0.5), tapers = tapers)
    @test all(x -> all(isfinite.(x)), getestimate(result))

    # Test large radii
    large_radii = [2.0, 100.0]
    @test_throws ArgumentError centered_l_function(
        points_data, radii = large_radii, nk = (8, 8), kmax = (0.5, 0.5), tapers = tapers)

    # Test centered L function at zero (should be zero or near zero)
    zero_radius = [0.0]
    result_zero = centered_l_function(
        points_data, radii = zero_radius, nk = (8, 8), kmax = (0.5, 0.5), tapers = tapers)
    @test abs(getestimate(result_zero)[1]) < 1e-6  # centered_L(0) should be approximately 0
end

@testset "in place tests $return_type" for return_type in [:single, :tuple, :vector]
    points_data = make_points_example(rng, n_processes = 1, dim = 2,
        return_type = :single, point_number = 20)
    tapers = sin_taper_family((3, 3), getregion(points_data))
    radii = [0.1, 0.3, 0.5]
    spectrum = spectra(points_data, nk = (16, 16), kmax = (0.5, 0.5), tapers = tapers)
    c_fun = c_function(spectrum, radii = radii)
    k_fun = k_function(spectrum, radii = radii)
    l_fun = l_function(spectrum, radii = radii)

    # from spatial data
    alloc = @ballocated centered_l_function!($spectrum, radii = $radii) samples=1
    @test_broken alloc == 0

    # from c function
    alloc = @ballocated centered_l_function!($c_fun) samples=1
    @test alloc == 0

    # from k function
    alloc = @ballocated centered_l_function!($k_fun) samples=1
    @test alloc == 0

    # from l function
    alloc = @ballocated centered_l_function!($l_fun) samples=1
    @test alloc == 0
end
