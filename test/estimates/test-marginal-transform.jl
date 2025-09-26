using SpatialMultitaper, Test, StableRNGs, StaticArrays
include("../test_utilities/TestUtils.jl")
using .TestUtils

import SpatialMultitaper: MarginallyTransformedEstimate, apply_marginal_transform,
                          getestimatename, getoriginaltype, gettransformname, getestimate,
                          getargument, Spectra, getprocessinformation,
                          getestimationinformation

@testset "MarginallyTransformedEstimate Construction" begin
    rng = StableRNG(123)

    @testset "Basic construction" begin
        # Create base spectra
        data, region = make_points_example(rng, n_processes = 2, point_number = 40)
        spec = spectra(data, region, nfreq = (6, 6), fmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        # Apply transform
        real_spec = real(spec)

        @test real_spec isa MarginallyTransformedEstimate
        @test getargument(real_spec) == getargument(spec)
        @test getprocessinformation(real_spec) == getprocessinformation(spec)
        @test getestimationinformation(real_spec) == getestimationinformation(spec)
        @test size(real_spec) == size(spec)
        @test embeddim(real_spec) == embeddim(spec)
    end
end

@testset "Mathematical Transformations" begin
    rng = StableRNG(123)
    data, region = make_points_example(rng, n_processes = 2, point_number = 30)
    spec = spectra(data, region, nfreq = (4, 4), fmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    @testset "real transform" begin
        real_spec = real(spec)

        # Check that imaginary parts are removed
        estimate = getestimate(real_spec)
        original_estimate = getestimate(spec)

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
        estimate = getestimate(imag_spec)
        original_estimate = getestimate(spec)

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
        estimate = getestimate(abs_spec)
        original_estimate = getestimate(spec)

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
        estimate = getestimate(abs2_spec)

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
        estimate = getestimate(conj_spec)
        original_estimate = getestimate(spec)

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
        estimate = getestimate(log_spec)
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
        estimate = getestimate(exp_spec)
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

        result = apply_marginal_transform(real, data)
        @test size(result) == size(data)
        @test all(x -> all(imag.(x) .≈ 0), result)
    end

    @testset "Regular array input" begin
        data = rand(ComplexF64, 3, 3, 5, 5)  # P x Q x freq1 x freq2

        result = apply_marginal_transform(abs, data)
        @test size(result) == size(data)
        @test all(real.(result) .≥ 0)
        @test all(imag.(result) .≈ 0)
    end
end

@testset "Name and Type Information" begin
    rng = StableRNG(123)
    data, region = make_points_example(rng, n_processes = 2, point_number = 30)
    spec = spectra(data, region, nfreq = (4, 4), fmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    @testset "Estimate naming" begin
        real_spec = real(spec)
        name = getestimatename(real_spec)

        # Should combine transform name with original estimate name
        @test occursin("real", lowercase(name))
        @test occursin("spectra", lowercase(name))  # Original type name
    end

    @testset "Type information" begin
        abs_spec = abs(spec)
        @test getoriginaltype(typeof(abs_spec)) <: Spectra
        @test gettransformname(typeof(abs_spec)) == "abs"  # Stripped function name
    end
end

@testset "Chaining Transformations" begin
    rng = StableRNG(123)
    data, region = make_points_example(rng, n_processes = 2, point_number = 30)
    spec = spectra(data, region, nfreq = (4, 4), fmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    @testset "Multiple transforms" begin
        # Chain: spec -> abs -> log
        abs_spec = abs(spec)
        log_abs_spec = log(abs_spec)

        @test log_abs_spec isa MarginallyTransformedEstimate
        @test size(log_abs_spec) == size(spec)

        # Check that the transformation actually happened
        estimate = getestimate(log_abs_spec)
        if estimate isa AbstractArray{<:SMatrix}
            @test all(x -> all(isreal(x) for x in x), estimate)
        else
            @test all(isreal.(estimate))
        end
    end
end

@testset "Indexing Transformed Estimates" begin
    rng = StableRNG(123)
    data, region = make_points_example(rng, n_processes = 3, point_number = 30)
    spec = spectra(data, region, nfreq = (4, 4), fmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    real_spec = real(spec)

    @testset "Process indexing" begin
        real_sub = real_spec[1, 2]
        @test real_sub isa MarginallyTransformedEstimate
        @test size(real_sub) == (1, 1)
        @test getargument(real_sub) == getargument(real_spec)
    end

    @testset "Frequency indexing" begin
        real_freq = real_spec[1, 1, 2, 3]
        @test real_freq isa MarginallyTransformedEstimate
        @test size(real_freq) == (1, 1)

        # Check frequency was correctly indexed
        freq_arg = getargument(real_freq)
        @test length(freq_arg[1]) == 1
        @test length(freq_arg[2]) == 1
    end
end

@testset "Type Stability and Performance" begin
    rng = StableRNG(123)
    data, region = make_points_example(rng, n_processes = 2, point_number = 20)
    spec = spectra(data, region, nfreq = (4, 4), fmax = (0.2, 0.2),
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
        real_est = getestimate(real_spec)
        if real_est isa AbstractArray{<:SMatrix}
            @test eltype(eltype(real_est)) <: Real
        else
            @test eltype(real_est) <: Real
        end

        # Abs transform should also produce real output
        abs_est = getestimate(abs_spec)
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
        data, region = make_points_example(rng, n_processes = 1, point_number = 10)
        spec = spectra(data, region, nfreq = (2, 2), fmax = (0.1, 0.1),
            tapers = sin_taper_family((1, 1), region))

        abs_spec = abs(spec)
        @test size(abs_spec) == (1, 1)
        @test embeddim(abs_spec) == 2
    end
end
