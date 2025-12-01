using SpatialMultitaper, Test, StableRNGs, StaticArrays
include("../../test_utilities/TestUtils.jl")
using .TestUtils

import SpatialMultitaper: MarginallyTransformedEstimate, apply_marginal_transform!,
                          get_estimate_name, get_original_type, get_transform_name,
                          get_estimates,
                          get_evaluation_points, Spectra, get_process_information,
                          get_estimation_information

@testset "MarginallyTransformedEstimate Construction" begin
    rng = StableRNG(123)

    @testset "Basic construction" begin
        # Create base spectra
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 40)
        region = getregion(data)
        spec = spectra(data, nk = (6, 6), kmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        # Apply transform
        real_spec = real(spec)

        @test real_spec isa MarginallyTransformedEstimate
        @test get_evaluation_points(real_spec) == get_evaluation_points(spec)
        @test get_process_information(real_spec) == get_process_information(spec)
        @test get_estimation_information(real_spec) == get_estimation_information(spec)
        @test size(real_spec) == size(spec)
        @test embeddim(real_spec) == embeddim(spec)
    end
end

@testset "Mathematical Transformations" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)
    region = getregion(data)
    spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    @testset "real transform" begin
        real_spec = real(spec)

        # Check that imaginary parts are removed
        estimate = get_estimates(real_spec)
        original_estimate = get_estimates(spec)

        if original_estimate isa AbstractArray{<:SMatrix}
            # For SMatrix case, check element-wise
            @test all(x -> all(imag.(x) .≈ 0), estimate)
            @test all(real.(x) ≈ y for (x, y) in zip(original_estimate, estimate))
        else
            # For regular array case
            @test all(imag.(estimate) .≈ 0)
            @test real.(original_estimate) ≈ estimate
        end
    end

    @testset "imag transform" begin
        imag_spec = imag(spec)
        estimate = get_estimates(imag_spec)
        original_estimate = get_estimates(spec)

        if original_estimate isa AbstractArray{<:SMatrix}
            @test all(x -> all(imag.(x) .≈ 0), estimate)  # Result should be real
            @test all(imag.(y) ≈ x for (x, y) in zip(estimate, original_estimate))
        else
            @test all(imag.(estimate) .≈ 0)  # Result should be real
            @test imag.(original_estimate) ≈ estimate
        end
    end

    @testset "abs transform" begin
        abs_spec = abs(spec)
        estimate = get_estimates(abs_spec)
        original_estimate = get_estimates(spec)

        # All values should be non-negative and real
        if original_estimate isa AbstractArray{<:SMatrix}
            @test all(x -> all(real(x) ≥ 0 for x in x), estimate)
            @test all(x -> all(imag(x) ≈ 0 for x in x), estimate)
        else
            @test all(real.(estimate) .≥ 0)
            @test all(imag.(estimate) .≈ 0)
        end
    end

    @testset "abs2 transform" begin
        abs2_spec = abs2(spec)
        estimate = get_estimates(abs2_spec)

        # Should be real and non-negative
        if estimate isa AbstractArray{<:SMatrix}
            @test all(x -> all(real.(x) .≥ 0), estimate)
            @test all(x -> all(imag.(x) .≈ 0), estimate)
        else
            @test all(real.(estimate) .≥ 0)
            @test all(imag.(estimate) .≈ 0)
        end
    end

    @testset "conj transform" begin
        conj_spec = conj(spec)
        estimate = get_estimates(conj_spec)
        original_estimate = get_estimates(spec)

        if original_estimate isa AbstractArray{<:SMatrix}
            @test all(conj.(y) ≈ x for (x, y) in zip(estimate, original_estimate))
        else
            @test conj.(original_estimate) ≈ estimate
        end
    end

    @testset "log transform" begin
        # Use abs to ensure positive values for log
        abs_spec = abs(spec)
        log_spec = log(abs_spec)

        # Should be real (since input was positive)
        estimate = get_estimates(log_spec)
        if estimate isa AbstractArray{<:SMatrix}
            @test all(x -> all(isreal(x) for x in x), estimate)
        else
            @test all(isreal.(estimate))
        end
    end

    @testset "exp transform" begin
        real_spec = real(spec)  # Use real part to avoid complex exponentials
        exp_spec = exp(real_spec)

        # Exponential should always be positive
        estimate = get_estimates(exp_spec)
        if estimate isa AbstractArray{<:SMatrix}
            @test all(x -> all(real(x) > 0 for x in x), estimate)
        else
            @test all(real.(estimate) .> 0)
        end
    end
end

@testset "apply_marginal_transform function" begin
    rng = StableRNG(123)

    @testset "SMatrix input" begin
        # Create array of SMatrix
        using StaticArrays
        data = [SMatrix{2, 2}(rand(ComplexF64, 2, 2)) for _ in 1:5, _ in 1:5]
        out = similar(data)
        result = apply_marginal_transform!(conj, out, data)
        @test size(result) == size(data)
    end

    @testset "Regular array input" begin
        data = rand(ComplexF64, 3, 3, 5, 5)  # P x Q x k1 x k2
        out = zeros(Float64, 3, 3, 5, 5)
        result = apply_marginal_transform!(abs, out, data)
        @test size(result) == size(data)
        @test all(real.(result) .≥ 0)
        @test all(imag.(result) .≈ 0)
    end
end

@testset "Name and Type Information" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)
    region = getregion(data)
    spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    @testset "Estimate naming" begin
        real_spec = real(spec)
        name = get_estimate_name(real_spec)

        # Should combine transform name with original estimate name
        @test occursin("real", lowercase(name))
        @test occursin("spectra", lowercase(name))  # Original type name
    end

    @testset "Type information" begin
        abs_spec = abs(spec)
        @test get_original_type(typeof(abs_spec)) <: Spectra
        @test get_transform_name(abs_spec) == "abs"  # Stripped function name
    end
end

@testset "Chaining Transformations" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)
    region = getregion(data)
    spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    @testset "Multiple transforms" begin
        # Chain: spec -> abs -> log
        abs_spec = abs(spec)
        log_abs_spec = log(abs_spec)

        @test log_abs_spec isa MarginallyTransformedEstimate
        @test size(log_abs_spec) == size(spec)

        # Check that the transformation actually happened
        estimate = get_estimates(log_abs_spec)
        if estimate isa AbstractArray{<:SMatrix}
            @test all(x -> all(isreal(x) for x in x), estimate)
        else
            @test all(isreal.(estimate))
        end
    end
end

@testset "Indexing Transformed Estimates" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 3, return_type = :tuple, point_number = 30)
    region = getregion(data)
    spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    real_spec = real(spec)

    @testset "Process indexing" begin
        real_sub = real_spec[1, 2]
        @test real_sub isa MarginallyTransformedEstimate
        @test size(real_sub) == ()
        @test get_evaluation_points(real_sub) == get_evaluation_points(real_spec)
    end

    @testset "Wavenumber indexing" begin
        real_wavenumber = real_spec[1, 1, 2, 3]
        @test real_wavenumber isa MarginallyTransformedEstimate
        @test size(real_wavenumber) == ()

        # Check wavenumber was correctly indexed
        wavenumber_arg = get_evaluation_points(real_wavenumber)
        @test length(wavenumber_arg[1]) == 1
        @test length(wavenumber_arg[2]) == 1
    end
end

@testset "Type Stability and Performance" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 20)
    region = getregion(data)
    spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    @testset "Type stability" begin
        @test typeof(real(spec)) <: MarginallyTransformedEstimate
        @test typeof(abs(spec)) <: MarginallyTransformedEstimate
        @test typeof(log(abs(spec))) <: MarginallyTransformedEstimate
    end

    @testset "Preserve element types appropriately" begin
        complex_spec = spec  # Already complex
        real_spec = real(spec)
        abs_spec = abs(spec)

        # Real transform should produce real output
        real_est = get_estimates(real_spec)
        if real_est isa AbstractArray{<:SMatrix}
            @test eltype(eltype(real_est)) <: Real
        else
            @test eltype(real_est) <: Real
        end

        # Abs transform should also produce real output
        abs_est = get_estimates(abs_spec)
        if abs_est isa AbstractArray{<:SMatrix}
            @test eltype(eltype(abs_est)) <: Real
        else
            @test eltype(abs_est) <: Real
        end
    end
end

@testset "Edge Cases" begin
    rng = StableRNG(123)

    @testset "Single element transform" begin
        data = make_points_example(
            rng, n_processes = 1, return_type = :single, point_number = 10)
        region = getregion(data)
        spec = spectra(data, nk = (2, 2), kmax = (0.1, 0.1),
            tapers = sin_taper_family((1, 1), region))

        abs_spec = abs(spec)
        @test size(abs_spec) == ()
        @test embeddim(abs_spec) == 2
    end
end
