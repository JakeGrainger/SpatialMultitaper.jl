using Test, SpatialMultitaper, StableRNGs, StaticArrays, LinearAlgebra, QuadGK,
      SpecialFunctions
include("../../test_utilities/TestData.jl")
using .TestData

import SpatialMultitaper: get_estimates, get_evaluation_points, CFunction,
                          default_rotational_kernel,
                          getregion, _isotropic_c_weight, _anisotropic_c_weight,
                          _isotropic_c_weight_generic, _anisotropic_c_weight_generic,
                          IsotropicEstimate

#

rng = StableRNG(123)

@testset "basic pipeline tests" begin
    @testset "single" begin
        data = make_points_example(
            rng, n_processes = 1, dim = 2, point_number = 50, return_type = :single)
        tapers = sin_taper_family((3, 3), getregion(data))
        c = c_function(data; radii = 0:0.1:1, kmax = 2.0, tapers = tapers)
        @test get_estimates(c) isa AbstractVector{<:Real}
        @test c isa IsotropicEstimate
        @test abs(c) isa IsotropicEstimate
    end
    @testset "vector" begin
        data = make_points_example(
            rng, n_processes = 1, dim = 2, point_number = 50, return_type = :vector)
        tapers = sin_taper_family((3, 3), getregion(data))
        c = c_function(data; radii = 0:0.1:1, kmax = 2.0, tapers = tapers)
        @test get_estimates(c) isa AbstractArray{<:Real, 3}
        @test c isa IsotropicEstimate
        @test abs(c) isa IsotropicEstimate
    end
    @testset "tuple" begin
        data = make_points_example(
            rng, n_processes = 1, dim = 2, point_number = 50, return_type = :tuple)
        tapers = sin_taper_family((3, 3), getregion(data))
        c = c_function(data; radii = 0:0.1:1, kmax = 2.0, tapers = tapers)
        @test get_estimates(c) isa AbstractVector{<:SMatrix{1, 1, <:Real, 1}}
        @test c isa IsotropicEstimate
        @test abs(c) isa IsotropicEstimate
    end
end

# loop over 1d, 2d, 3d
@testset "Dimension $dim tests" for dim in 1:3

    # - c function from raw data is same as from Spatial data
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

        raw_result = c_function(
            points_raw, region, radii = radii, nk = nk, kmax = kmax, tapers = tapers)
        spatial_result = c_function(
            points_spatial, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

        @test get_evaluation_points(raw_result) ≈ get_evaluation_points(spatial_result)
        @test get_estimates(raw_result) ≈ get_estimates(spatial_result)
    end

    # - loop over single, tuple and vector
    @testset "Return type: $return_type" for return_type in [:single, :tuple, :vector]
        n_processes = return_type == :single ? 1 : 3

        # - - c function from SpatialData
        @testset "C-function from SpatialData" begin
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

            result = c_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test result isa CFunction
            @test get_evaluation_points(result) ≈ radii
            if return_type == :vector
                @test size(get_estimates(result)) ==
                      (n_processes, n_processes, length(radii))
            else
                @test length(get_estimates(result)) == length(radii)
            end
            @test all(x -> all(isfinite.(x)), get_estimates(result))
        end

        # - - c function from spectra
        @testset "C-function from spectra" begin
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

            # Then compute c_function from spectra
            c_from_spectra = c_function(spectrum, radii = radii)

            # Compare with direct computation
            c_direct = c_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test get_evaluation_points(c_from_spectra) ≈ get_evaluation_points(c_direct)
            @test get_estimates(c_from_spectra) ≈ get_estimates(c_direct)
        end

        # - - partial c function from SpatialData
        @testset "Partial c-function from SpatialData" begin
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

                partial_result = partial_c_function(points_data, radii = radii,
                    nk = nk, kmax = kmax, tapers = tapers)

                @test partial_result isa CFunction
                @test get_evaluation_points(partial_result) ≈ radii
                if return_type == :vector
                    @test size(get_estimates(partial_result)) ==
                          (n_processes, n_processes, length(radii))
                else
                    @test length(get_estimates(partial_result)) == length(radii)
                end
            end
        end

        # - - partial c function from partial spectra
        @testset "Partial c-function from partial spectra" begin
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

                # Then c_function from partial spectra
                c_from_partial = c_function(partial_spec, radii = radii)

                @test c_from_partial isa CFunction
                @test get_evaluation_points(c_from_partial) ≈ radii
            end
        end

        # - - partial c function from spectra
        @testset "Partial c-function from spectra" begin
            if n_processes > 1
                points_data = make_points_example(
                    rng, n_processes = n_processes, dim = dim,
                    return_type = return_type, point_number = 20)

                radii = dim == 1 ? [0.1] : [0.15]
                nk = dim == 1 ? (10,) : dim == 2 ? (6, 6) : (4, 4, 4)
                kmax = dim == 1 ? (0.6,) : dim == 2 ? (0.5, 0.5) : (0.3, 0.3, 0.3)
                bw = ntuple(_ -> 3, dim)
                region = getregion(points_data)
                tapers = sin_taper_family(bw, region)

                # Compute regular spectra first
                spectrum = spectra(
                    points_data, nk = nk, kmax = kmax, tapers = tapers)

                # Then partial c_function from regular spectra
                partial_from_spec = partial_c_function(spectrum, radii = radii)

                @test partial_from_spec isa CFunction
                @test get_evaluation_points(partial_from_spec) ≈ radii
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

            result = c_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            @test result isa CFunction
            @test get_evaluation_points(result) isa AbstractVector
            @test length(get_evaluation_points(result)) == length(radii)
            if return_type == :vector
                @test get_estimates(result) isa AbstractArray
                @test size(get_estimates(result)) ==
                      (n_processes, n_processes, length(radii))
            elseif return_type == :tuple
                @test get_estimates(result) isa AbstractVector
                @test length(get_estimates(result)) == length(radii)
                @test eltype(get_estimates(result)) <: SMatrix{n_processes, n_processes}
            else # return_type == :single
                @test get_estimates(result) isa AbstractVector
                @test length(get_estimates(result)) == length(radii)
            end

            # Check element types
            @test eltype(get_evaluation_points(result)) <: Number
            if return_type == :tuple
                @test eltype(get_estimates(result)) <: SMatrix
            else
                @test eltype(get_estimates(result)) <: Number
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

            result = c_function(
                points_data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)

            # Test indexing into results
            radii_result = get_evaluation_points(result)
            values_result = get_estimates(result)

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
    end
end

@testset "No radii given" begin
    points_data = make_points_example(rng, n_processes = 1, dim = 2,
        return_type = :single, point_number = 20)
    @test c_function(points_data, kmax = 0.5) isa CFunction
end

@testset "Edge cases and robustness" begin
    # Test very small radii
    points_data = make_points_example(rng, n_processes = 1, dim = 2,
        return_type = :single, point_number = 20)
    small_radii = [0.01, 0.05]
    bw = (3, 3)
    region = getregion(points_data)
    tapers = sin_taper_family(bw, region)
    result = c_function(
        points_data, radii = small_radii, nk = (8, 8), kmax = (0.5, 0.5), tapers = tapers)
    @test all(x -> all(isfinite.(x)), get_estimates(result))

    # radii outside region size should error
    large_radii = [2.0, 100.0]
    @test_throws ArgumentError c_function(
        points_data, radii = large_radii, nk = (8, 8), kmax = (0.5, 0.5), tapers = tapers)
end

@testset "isotropic transformation in 2d" begin
    dim = 2

    # Define test parameters
    radii = dim == 1 ? [0.1, 0.5, 1.0] :
            dim == 2 ? [0.1, 0.3, 0.8] : [0.1, 0.2, 0.5]
    nk = dim == 1 ? (32,) : dim == 2 ? (16, 16) : (8, 8, 8)
    kmax = dim == 1 ? (2.0,) : dim == 2 ? (1.5, 1.5) : (1.0, 1.0, 1.0)
    bw = ntuple(_ -> 3, dim)
    points_raw, region = make_points_example(rng, n_processes = 1, dim = dim,
        return_type = :raw, point_number = 50)
    points_single = spatial_data(points_raw[1], region)
    points_vec = spatial_data([points_raw[1]], region)
    points_tuple = spatial_data((points_raw[1],), region)

    tapers = sin_taper_family(bw, region)
    c_single = c_function(
        points_single, radii = radii, nk = nk, kmax = kmax, tapers = tapers,
        rotational_method = default_rotational_kernel(nk, kmax))
    c_vec = c_function(
        points_vec, radii = radii, nk = nk, kmax = kmax, tapers = tapers,
        rotational_method = default_rotational_kernel(nk, kmax))
    c_tuple = c_function(
        points_tuple, radii = radii, nk = nk, kmax = kmax, tapers = tapers,
        rotational_method = default_rotational_kernel(nk, kmax))

    @test c_single isa CFunction
    @test c_vec isa CFunction
    @test c_tuple isa CFunction

    @test get_estimates(c_single) ≈ get_estimates(c_vec)[1, 1, :]
    @test get_estimates(c_single) ≈ getindex.(get_estimates(c_tuple), 1, 1)
end

R = [0.0, 0.1, 1.0, 10.0]
U = [0.1, 0.5, 1.0, 10.0] # with the SL only makes sense for u >= sl/2
SL = [0.1]

@testset "anisotropic transforms" begin
    @testset "r=$r, ||u||=$u, dim=$dim" for r in R, u in [0.0; U], dim in 1:3
        theoretical = _anisotropic_c_weight(r, u, Val{dim}())
        general = _anisotropic_c_weight_generic(r, u, Val{dim}())
        @test theoretical≈general atol=1e-6
    end
end
@testset "anisotropic integrals" begin
    @testset "r=$r, ||u||=$u, dim=$dim, sl=$sl" for r in R, u in U, sl in SL, dim in 1:3
        theoretical = _isotropic_c_weight(r, u, sl, Val{dim}())
        general = _isotropic_c_weight_generic(r, u, sl, Val{dim}())
        if r == u == 10
            @test_broken theoretical≈general atol=1e-5
        else
            @test theoretical≈general atol=1e-5
        end
    end
end

@testset "1d integration tests" begin
    @testset "r=$r, ||u||=$u, dim=$dim, sl=$sl" for r in R, u in U, sl in SL, dim in 1:3
        numerical = 2 * pi^(dim / 2) / gamma(dim / 2) *
                    quadgk(x -> x^(dim - 1) * _anisotropic_c_weight(r, x, Val{dim}()),
            u - sl / 2, u + sl / 2)[1]

        theoretical = _isotropic_c_weight(r, u, sl, Val{dim}())
        @test theoretical≈numerical atol=1e-5
    end
end

@testset "checking isotropic integrals" begin
    @testset "r=$r, ||u||=$u, dim=$dim, sl=$sl" for r in R, u in U, sl in SL, dim in 1:4
        # Value should be r^{d/2} * A_d * ∫_{k - sl/2}^{k + sl/2} x^{d/2 - 1} * J_{d/2}(2π r x) dx,
        # where A_d = 2π^{d/2} / Γ(d/2)
        theoretical = _isotropic_c_weight_generic(r, norm(u), sl, Val{dim}())
        int_approx = quadgk(
            x -> (x^(dim / 2 - 1) * besselj(dim / 2, 2π * r * x)), u - sl / 2, u + sl / 2)[1]
        numerical = 2 * (pi * r)^(dim / 2) / gamma(dim / 2) * int_approx
        if r == u == 10
            @test_broken theoretical≈numerical atol=1e-5
        else
            @test theoretical≈numerical atol=1e-5
        end
    end
end
