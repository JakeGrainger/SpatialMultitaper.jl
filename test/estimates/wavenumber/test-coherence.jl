using SpatialMultitaper, Test, StableRNGs, LinearAlgebra
include("../../test_utilities/TestUtils.jl")
using .TestUtils

import SpatialMultitaper: Coherence, MarginallyTransformedEstimate, getestimate,
                          getargument, gettransformtype, getprocessinformation,
                          getestimationinformation

@testset "coherence matrix function" begin
    @testset "2x2 matrix" begin
        # Create a known spectral matrix
        S = [4.0 2.0+1.0im; 2.0-1.0im 9.0]
        C = coherence(S)

        # Coherence should normalize by sqrt of diagonal elements
        @test C[1, 1] ≈ 1.0  # Auto-coherence is always 1
        @test C[2, 2] ≈ 1.0
        @test C[1, 2] ≈ (2.0 + 1.0im) / (sqrt(4.0) * sqrt(9.0))  # = (2+i)/6
        @test C[2, 1] ≈ conj(C[1, 2])  # Hermitian symmetry
    end

    @testset "3x3 matrix" begin
        S = [1.0 0.5im -0.3; -0.5im 2.0 0.7+0.2im; -0.3 0.7-0.2im 3.0]
        C = coherence(S)

        # Check diagonal elements
        @test all(isapprox(C[i, i], 1.0) for i in 1:3)

        # Check normalization
        @test C[1, 2] ≈ 0.5im / sqrt(2.0)
        @test C[1, 3] ≈ -0.3 / sqrt(3.0)
        @test C[2, 3] ≈ (0.7 + 0.2im) / sqrt(6.0)

        # Check Hermitian property
        @test C ≈ C'
    end
end

@testset "Coherence struct and construction" begin
    rng = StableRNG(123)

    @testset "From Spectra" begin
        # Create test spectra
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 40)
        nk = (6, 6)
        kmax = (0.3, 0.3)
        region = getregion(data)
        tapers = sin_taper_family((2, 2), region)

        spec = spectra(data, nk = nk, kmax = kmax, tapers = tapers)
        coh = coherence(spec)

        @test coh isa Coherence
        @test getargument(coh) == getargument(spec)  # Same frequencies
        @test size(coh) == size(spec)
        @test embeddim(coh) == embeddim(spec)
        @test getprocessinformation(coh) == getprocessinformation(spec)
        @test getestimationinformation(coh) == getestimationinformation(spec)
    end

    @testset "Direct from data" begin
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 40)
        nk = (6, 6)
        kmax = (0.3, 0.3)
        region = getregion(data)
        tapers = sin_taper_family((2, 2), region)

        coh = coherence(data, nk = nk, kmax = kmax, tapers = tapers)
        @test coh isa Coherence
        @test size(coh) == (2, 2)
    end
end

@testset "partial_coherence" begin
    @testset "Matrix function" begin
        # Test with known matrix
        S = [2.0 0.8; 0.8 1.0]
        S_inv = inv(S)
        C_partial = partial_coherence(S)
        C_expected = -coherence(S_inv)

        @test C_partial ≈ C_expected

        # For 2x2 case, check specific formula
        expected_12 = -S_inv[1, 2] / (sqrt(S_inv[1, 1]) * sqrt(S_inv[2, 2]))
        @test C_partial[1, 2] ≈ expected_12
    end

    @testset "From Spectra - Marginal" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 3, return_type = :tuple, point_number = 30)
        nk = (4, 4)
        kmax = (0.2, 0.2)
        region = getregion(data)
        tapers = sin_taper_family((2, 2), region)

        spec = spectra(data, nk = nk, kmax = kmax, tapers = tapers)
        coh_partial = partial_coherence(spec)

        @test coh_partial isa Coherence
        @test size(coh_partial) == size(spec)
    end

    @testset "From Spectra - Partial trait" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 30)
        region = getregion(data)

        # Create a partial spectrum first
        spec_partial = partial_spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        coh_from_partial = partial_coherence(spec_partial)

        # For partial spectra, partial_coherence should equal coherence
        coh_regular = coherence(spec_partial)
        @test getestimate(coh_from_partial) ≈ getestimate(coh_regular)
    end
end

@testset "magnitude_coherence" begin
    @testset "From Spectra" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 40)
        nk = (6, 6)
        kmax = (0.3, 0.3)
        region = getregion(data)
        tapers = sin_taper_family((2, 2), region)

        spec = spectra(data, nk = nk, kmax = kmax, tapers = tapers)
        mag_coh = magnitude_coherence(spec)

        # Should be real and positive
        estimate = getestimate(mag_coh)
        @test all(x -> all(real(x) ≥ 0 for x in x), estimate)  # All elements ≥ 0
        @test all(x -> all(imag(x) ≈ 0 for x in x), estimate)  # All elements real
    end

    @testset "From Coherence" begin
        rng = StableRNG(123)
        data = make_points_example(rng, n_processes = 2, return_type = :tuple)
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        coh = coherence(spec)
        mag_coh = magnitude_coherence(coh)

        # Test that it's the absolute value
        expected = map(x -> abs.(x), getestimate(coh))
        @test getestimate(mag_coh) ≈ expected
    end

    @testset "Direct from data" begin
        rng = StableRNG(123)
        data = make_points_example(rng, n_processes = 2, return_type = :tuple)
        region = getregion(data)
        mag_coh = magnitude_coherence(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        @test mag_coh isa MarginallyTransformedEstimate
        @test gettransformtype(mag_coh) === typeof(abs)
    end
end

@testset "magnitude_squared_coherence" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)
    region = getregion(data)
    spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    @testset "From Spectra" begin
        mag2_coh = magnitude_squared_coherence(spec)
        estimate = getestimate(mag2_coh)
        @test all(x -> all(real(x) ≥ 0 for x in x), estimate)  # All elements ≥ 0
        @test all(x -> all(imag(x) ≈ 0 for x in x), estimate)  # All elements real
    end

    @testset "From Coherence" begin
        coh = coherence(spec)
        mag2_coh = magnitude_squared_coherence(coh)
        expected = map(x -> abs2.(x), getestimate(coh))
        @test getestimate(mag2_coh) ≈ expected
    end
end

@testset "phase function" begin
    @testset "From Spectra" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 30)
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        phase_est = phase(spec)
        estimate = getestimate(phase_est)

        # Phase should be real and between -π and π
        @test all(x -> all(imag(x) ≈ 0 for x in x), estimate)  # All elements real
        @test all(x -> all(-π ≤ real(x) ≤ π for x in x), estimate)  # In [-π, π]
    end

    @testset "From Coherence" begin
        rng = StableRNG(123)
        data = make_points_example(rng, n_processes = 2, return_type = :tuple)
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        coh = coherence(spec)

        phase_from_coh = phase(coh)
        phase_from_spec = phase(spec)

        # Should give same result (both use angle)
        @test getestimate(phase_from_coh) ≈ getestimate(phase_from_spec)
    end

    @testset "Direct from data" begin
        rng = StableRNG(123)
        data = make_points_example(rng, n_processes = 2, return_type = :tuple)
        region = getregion(data)
        phase_direct = phase(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        @test phase_direct isa MarginallyTransformedEstimate
        @test gettransformtype(phase_direct) === typeof(angle)
    end
end

@testset "Mathematical Properties" begin
    rng = StableRNG(123)

    @testset "Coherence properties" begin
        # Create deterministic test matrix
        S = [26.0 2.0+im; 2.0-im 1.0]
        C = coherence(S)

        # Auto-coherence should be 1
        @test C[1, 1] ≈ 1.0
        @test C[2, 2] ≈ 1.0

        # Magnitude should be ≤ 1
        @test abs(C[1, 2]) ≤ 1.0
        @test abs(C[2, 1]) ≤ 1.0

        # Hermitian symmetry
        @test C[1, 2] ≈ conj(C[2, 1])
    end

    @testset "Relationship between transforms" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 50)
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        coh = coherence(spec)
        mag_coh = magnitude_coherence(spec)
        mag2_coh = magnitude_squared_coherence(spec)

        # |coherence|^2 should equal magnitude_squared_coherence
        @test map(x -> abs2.(x), getestimate(coh)) ≈ getestimate(mag2_coh)

        # |coherence| should equal magnitude_coherence
        @test map(x -> abs.(x), getestimate(coh)) ≈ getestimate(mag_coh)
    end
end

@testset "Edge Cases" begin
    @testset "Perfect coherence" begin
        # Identity matrix should give perfect coherence
        S = Matrix{ComplexF64}(I, 2, 2)
        C = coherence(S)

        @test C ≈ S  # Coherence of identity is identity
        @test abs(C[1, 2]) ≈ 0.0  # Off-diagonal should be zero
    end

    @testset "Single process" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 1, return_type = :single, point_number = 30)
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        coh = coherence(spec)
        @test size(coh) == (1, 1)

        # Single process coherence should be 1 everywhere
        estimate = getestimate(coh)
        @test all(x[1] ≈ 1.0 for x in estimate)
    end
end
