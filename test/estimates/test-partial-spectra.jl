using SpatialMultitaper, Test, StableRNGs, LinearAlgebra, StaticArrays
include("../test_utilities/TestUtils.jl")
using .TestUtils

import SpatialMultitaper: partial_spectra, partial_spectra_uncorrected, getestimate,
                          getargument, is_partial, MarginallyTransformedEstimate, Spectra,
                          getprocessinformation, getestimationinformation

@testset "partial_spectra matrix function" begin
    @testset "2x2 SMatrix" begin
        # Create a known spectral matrix
        S = SMatrix{2, 2, ComplexF64, 4}([2.0+0.0im 0.5+0.5im; 0.5-0.5im 1.0+0.0im])

        # Test uncorrected version first
        P_uncorrected = partial_spectra(S, nothing)

        # For 2x2 case, manual calculation
        S_inv = inv(S)
        expected_11 = 1 / S_inv[1, 1]
        expected_22 = 1 / S_inv[2, 2]
        expected_12 = -S_inv[1, 2] / (S_inv[1, 1] * S_inv[2, 2] - abs2(S_inv[1, 2]))
        expected_21 = conj(expected_12)

        @test P_uncorrected[1, 1] ≈ expected_11
        @test P_uncorrected[2, 2] ≈ expected_22
        @test P_uncorrected[1, 2] ≈ expected_12
        @test P_uncorrected[2, 1] ≈ expected_21
    end

    @testset "3x3 regular matrix" begin
        S = [3.0 1.0 0.5; 1.0 2.0 0.3; 0.5 0.3 1.5]
        P = partial_spectra(S, nothing)

        # Check dimensions
        @test size(P) == (3, 3)

        # Diagonal elements should be 1/g_ii where g_ii are diagonal of inv(S)
        S_inv = inv(S)
        for i in 1:3
            @test P[i, i] ≈ 1 / S_inv[i, i]
        end

        # Off-diagonal elements: -g_ij / (g_ii * g_jj - |g_ij|²)
        @test P[1, 2] ≈ -S_inv[1, 2] / (S_inv[1, 1] * S_inv[2, 2] - abs2(S_inv[1, 2]))
    end

    @testset "Bias correction with ntapers" begin
        S = SMatrix{2, 2, Float64, 4}([2.0 0.8; 0.8 1.0])
        ntapers = 10

        P_corrected = partial_spectra(S, ntapers)
        P_uncorrected = partial_spectra(S, nothing)

        # Check that correction factor is applied
        # Correction for off-diagonal: M/(M-Q+2) = 10/(10-2+2) = 1.0
        # Correction for diagonal: M/(M-Q+1) = 10/(10-2+1) = 10/9

        @test P_corrected[1, 1] ≈ P_uncorrected[1, 1] * 10 / 9
        @test P_corrected[2, 2] ≈ P_uncorrected[2, 2] * 10 / 9
        @test P_corrected[1, 2] ≈ P_uncorrected[1, 2] * 1.0
        @test P_corrected[2, 1] ≈ P_uncorrected[2, 1] * 1.0
    end
end

@testset "partial_spectra from Spectra" begin
    rng = StableRNG(123)

    @testset "Basic functionality" begin
        data = make_points_example(
            rng, n_processes = 3, return_type = :tuple, point_number = 40)
        nfreq = (6, 6)
        fmax = (0.3, 0.3)
        region = getregion(data)
        tapers = sin_taper_family((2, 2), region)

        # Create marginal spectra
        spec = spectra(data, nfreq = nfreq, fmax = fmax, tapers = tapers)

        # Convert to partial spectra
        partial_spec = partial_spectra(spec)

        @test partial_spec isa Spectra
        @test is_partial(partial_spec) == true  # Should have PartialTrait
        @test size(partial_spec) == size(spec)
        @test embeddim(partial_spec) == embeddim(spec)
        @test getargument(partial_spec) == getargument(spec)  # Same frequencies
        @test getprocessinformation(partial_spec) == getprocessinformation(spec)
        @test getestimationinformation(partial_spec) == getestimationinformation(spec)
    end

    @testset "Direct from data" begin
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 30)
        nfreq = (4, 4)
        fmax = (0.2, 0.2)
        region = getregion(data)
        tapers = sin_taper_family((2, 2), region)

        partial_spec = partial_spectra(data, nfreq = nfreq, fmax = fmax, tapers = tapers)

        @test partial_spec isa Spectra
        @test is_partial(partial_spec) == true
        @test size(partial_spec) == (2, 2)
    end

    @testset "Uncorrected version" begin
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 30)
        region = getregion(data)
        spec = spectra(data, nfreq = (4, 4), fmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        partial_corrected = partial_spectra(spec)
        partial_uncorrected = partial_spectra_uncorrected(spec)

        @test partial_uncorrected isa Spectra
        @test is_partial(partial_uncorrected) == true

        # Results should be different (uncorrected should be smaller in magnitude typically)
        @test getestimate(partial_corrected) ≢ getestimate(partial_uncorrected)
    end
end

@testset "Mathematical Properties" begin
    rng = StableRNG(123)

    @testset "Positive definite input" begin
        # Create a positive definite matrix
        A = [2.0 0.5; 0.5 1.0]
        @test isposdef(A)

        P = partial_spectra(A, nothing)

        # Diagonal elements should be positive (since A is positive definite)
        @test real(P[1, 1]) > 0
        @test real(P[2, 2]) > 0
    end

    @testset "Hermitian property" begin
        # Start with Hermitian matrix
        S = [2.0 1.0+0.5im; 1.0-0.5im 1.5]
        @test S ≈ S'

        P = partial_spectra(S, nothing)

        # Result should also be Hermitian
        @test P ≈ P'
    end

    @testset "Inverse relationship" begin
        # For 2x2 case, can check specific relationship with inverse
        S = [3.0 1.0; 1.0 2.0]  # Real symmetric positive definite
        P = partial_spectra(S, nothing)
        S_inv = inv(S)

        # Diagonal: P[i,i] = 1/S_inv[i,i]
        @test P[1, 1] ≈ 1 / S_inv[1, 1]
        @test P[2, 2] ≈ 1 / S_inv[2, 2]

        # Off-diagonal has more complex relationship involving denominators
        expected_12 = -S_inv[1, 2] / (S_inv[1, 1] * S_inv[2, 2] - S_inv[1, 2]^2)
        @test P[1, 2] ≈ expected_12
    end
end

@testset "Different Matrix Sizes" begin
    @testset "1x1 matrix (trivial case)" begin
        S = [2.0;;]  # 1x1 matrix
        P = partial_spectra(S, nothing)

        @test size(P) == (1, 1)
        @test P[1, 1] ≈ 2.0 # Should be same as input because nothing to condition on
    end

    @testset "4x4 matrix" begin
        rng = StableRNG(123)
        # Create a positive definite matrix
        A = rand(rng, 4, 4)
        S = A'A + I  # Ensure positive definite

        P = partial_spectra(S, nothing)

        @test size(P) == (4, 4)
        @test P ≈ P'  # Should be Hermitian

        # Check diagonal elements are positive
        @test all(real(P[i, i]) > 0 for i in 1:4)
    end
end

@testset "Bias Correction Analysis" begin
    @testset "Effect of ntapers parameter" begin
        S = SMatrix{3, 3, Float64, 9}([3.0 1.0 0.5; 1.0 2.0 0.3; 0.5 0.3 1.5])

        P_uncorrected = partial_spectra(S, nothing)
        P_corrected_5 = partial_spectra(S, 5)
        P_corrected_20 = partial_spectra(S, 20)

        # Higher ntapers should give results closer to uncorrected
        # (since correction factor approaches 1 as M increases)
        for i in 1:3, j in 1:3
            diff_5 = abs(P_corrected_5[i, j] - P_uncorrected[i, j])
            diff_20 = abs(P_corrected_20[i, j] - P_uncorrected[i, j])

            if i == j
                # Diagonal correction: M/(M-Q+1)
                @test P_corrected_5[i, j] ≈ P_uncorrected[i, j] * 5 / (5 - 3 + 1)
                @test P_corrected_20[i, j] ≈ P_uncorrected[i, j] * 20 / (20 - 3 + 1)
            else
                # Off-diagonal correction: M/(M-Q+2)
                @test P_corrected_5[i, j] ≈ P_uncorrected[i, j] * 5 / (5 - 3 + 2)
                @test P_corrected_20[i, j] ≈ P_uncorrected[i, j] * 20 / (20 - 3 + 2)
            end
        end
    end
end

@testset "Integration with Estimate Framework" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)
    region = getregion(data)
    spec = spectra(data, nfreq = (4, 4), fmax = (0.2, 0.2),
        tapers = sin_taper_family((2, 2), region))

    partial_spec = partial_spectra(spec)

    @testset "Indexing operations" begin
        # Test process indexing
        partial_sub = partial_spec[1, 2]
        @test size(partial_sub) == (1, 1)
        @test is_partial(partial_sub) == true

        # Test frequency indexing
        partial_freq = partial_spec[1, 1, 2, 3]
        @test size(partial_freq) == (1, 1)
        @test is_partial(partial_freq) == true
    end

    @testset "Compatibility with other transforms" begin
        # Should be able to apply other transforms to partial spectra
        abs_partial = abs(partial_spec)
        @test abs_partial isa MarginallyTransformedEstimate

        real_partial = real(partial_spec)
        @test real_partial isa MarginallyTransformedEstimate
    end
end

@testset "Numerical Stability" begin
    @testset "Near-singular matrix" begin
        # Create a matrix that's close to singular
        S = [1.0 0.999; 0.999 1.0]
        @test cond(S) > 100  # Ill-conditioned

        P = partial_spectra(S, nothing)

        # Should still produce finite results
        @test all(isfinite, P)
        @test P ≈ P'  # Should maintain Hermitian property
    end

    @testset "Very small values" begin
        S = [1e-6 1e-8; 1e-8 1e-6]
        P = partial_spectra(S, nothing)

        @test all(isfinite, P)
        @test all(!isnan, P)
    end

    @testset "Large values" begin
        S = [1e6 1e4; 1e4 1e6]
        P = partial_spectra(S, nothing)

        @test all(isfinite, P)
        @test all(!isnan, P)
    end
end

@testset "Edge Cases and Error Conditions" begin
    @testset "Zero matrix" begin
        S = zeros(2, 2)

        # Should throw error or handle gracefully (matrix is not invertible)
        @test_throws SingularException partial_spectra(S, nothing)
    end

    @testset "Identity matrix" begin
        S = Matrix{Float64}(I, 3, 3)
        P = partial_spectra(S, nothing)

        # Partial spectra of identity should be identity
        # Since inv(I) = I, diagonal elements are 1, off-diagonal are 0
        @test P ≈ S
    end

    @testset "Single-process spectra edge case" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 1, return_type = :single, point_number = 20)
        region = getregion(data)
        spec = spectra(data, nfreq = (4, 4), fmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        partial_spec = partial_spectra(spec)

        @test size(partial_spec) == (1, 1)
        @test is_partial(partial_spec) == true

        # For single process, partial spectrum should equal to the original
        estimate = getestimate(partial_spec)
        original_estimate = getestimate(spec)
        @test original_estimate ≈ estimate
    end
end
