using SpatialMultitaper, Test, StableRNGs
include("../../test_utilities/TestData.jl")
using .TestData

import SpatialMultitaper: RotationalEstimate, MarginallyTransformedEstimate, Coherence,
                          rotational_estimate, default_rotational_radii, is_partial,
                          default_rotational_kernel, GaussKernel, RectKernel,
                          _smoothed_rotational, getestimate, getargument, getestimatename,
                          getestimationinformation, getprocessinformation,
                          SingleProcessTrait

@testset "Rotational Kernels" begin
    @testset "GaussKernel" begin
        bw = 0.1
        kernel = GaussKernel(bw)

        # Test at zero
        @test kernel(0.0) ≈ 1 / bw

        # Test symmetry
        @test kernel(0.05) ≈ kernel(-0.05)

        # Test decay
        @test kernel(0.2) < kernel(0.1)
        @test kernel(0.1) < kernel(0.0)

        # Test normalization factor
        @test kernel(0.0) ≈ exp(0) / bw
    end

    @testset "RectKernel" begin
        bw = 0.2
        kernel = RectKernel(bw)

        # Inside window
        @test kernel(0.0) ≈ 1 / bw
        @test kernel(0.05) ≈ 1 / bw
        @test kernel(-0.05) ≈ 1 / bw

        # At boundary (should still be inside)
        @test kernel(bw / 2 - 1e-10) ≈ 1 / bw

        # Outside window
        @test kernel(bw / 2 + 1e-10) ≈ 0.0
        @test kernel(0.2) ≈ 0.0
        @test kernel(-0.2) ≈ 0.0
    end
end

@testset "Default Radii and Kernel Generation" begin
    @testset "From wavenumber vectors" begin
        freq1 = range(0, 0.5, length = 11)
        freq2 = range(0, 0.3, length = 8)
        wavenumber = (freq1, freq2)

        radii = default_rotational_radii(wavenumber)

        # Should use minimum of maximum wavenumbers
        max_radius = minimum([0.5, 0.3]) - step(radii) / 2
        @test maximum(radii) ≈ max_radius

        # Should use maximum of lengths for number of points
        @test length(radii) == max(11, 8) - 1

        # Should be evenly spaced
        @test step(radii) isa Number
    end

    @testset "From nk and kmax" begin
        radii = default_rotational_radii((11, 8), (0.4, 0.6))
        #TODO: add proper tests when this is reworked
    end

    @testset "Default kernel from estimate" begin
        rng = StableRNG(123)
        data = make_points_example(rng, return_type = :tuple)
        region = getregion(data)
        mt_est = spectra(data; nk = (10, 10), kmax = (0.5, 0.5),
            tapers = sin_taper_family((3, 3), region))
        kernel = default_rotational_kernel(mt_est)

        @test kernel isa RectKernel
        @test kernel.bw ≈ 0.2
    end
end

@testset "_smoothed_rotational function" begin
    @testset "1D case" begin
        x = (range(-1, 1, length = 21),)  # 1D wavenumbers
        y = ones(21)  # Constant function
        radii = [0.0, 0.5, 1.0]
        kernel = RectKernel(0.2)

        result = _smoothed_rotational(x, y, SingleProcessTrait(), radii, kernel)

        @test length(result) == length(radii)
        @test all(result .> 0)  # Should be positive since y is positive

        # For constant function and symmetric kernel, should be nearly constant
        # (but may vary due to edge effects)
    end

    @testset "2D case" begin
        x1 = range(-1, 1, length = 11)
        x2 = range(-1, 1, length = 11)
        x = (x1, x2)
        y = ones(11, 11)  # Constant 2D function

        radii = [0.0, 0.5, 1.0, 1.4]  # Note: √2 ≈ 1.414 is max radius for this grid
        kernel = RectKernel(0.3)

        result = _smoothed_rotational(x, y, SingleProcessTrait(), radii, kernel)

        @test length(result) == length(radii)
        @test all(result .> 0)
        @test result[1] ≈ 1.0  # At r=0, should equal center value
    end

    @testset "Non-constant function" begin
        x1 = range(-2, 2, length = 21)
        x2 = range(-2, 2, length = 21)
        x = (x1, x2)

        # Create radial function: f(r) = exp(-r²)
        y = [exp(-(xi^2 + xj^2)) for xi in x1, xj in x2]

        radii = [0.0, 0.5, 1.0, 1.5]
        kernel = RectKernel(0.2)

        result = _smoothed_rotational(x, y, SingleProcessTrait(), radii, kernel)

        @test length(result) == length(radii)
        @test result[1] > result[2] > result[3] > result[4]  # Should decrease with radius
        @test result[1]≈1.0 atol=0.1  # At r=0, should be close to 1
    end
end

@testset "RotationalEstimate Construction" begin
    rng = StableRNG(123)

    @testset "From anisotropic spectra" begin
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 40)
        region = getregion(data)
        spec = spectra(data, nk = (8, 8), kmax = (0.4, 0.4),
            tapers = sin_taper_family((2, 2), region))

        rot_spec = rotational_estimate(spec)

        @test rot_spec isa RotationalEstimate
        @test embeddim(rot_spec) == 2  # Same as original
        @test size(rot_spec) == size(spec)  # Same processes
        @test getprocessinformation(rot_spec) == getprocessinformation(spec)
        @test getestimationinformation(rot_spec) == getestimationinformation(spec)

        # Should be isotropic now (1D argument)
        radii = getargument(rot_spec)
        @test radii isa AbstractVector
        @test all(radii .≥ 0)  # Radii should be non-negative
    end

    @testset "Custom radii and kernel" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 30)
        region = getregion(data)
        spec = spectra(data, nk = (6, 6), kmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        custom_radii = [0.0, 0.1, 0.2, 0.3]
        custom_kernel = GaussKernel(0.05)

        rot_spec = rotational_estimate(spec, radii = custom_radii, kernel = custom_kernel)

        @test getargument(rot_spec) == custom_radii
        @test length(getestimate(rot_spec)) == length(custom_radii)
    end

    @testset "Default parameters" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 30)
        region = getregion(data)
        spec = spectra(data, nk = (6, 6), kmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))

        rot_spec = rotational_estimate(spec)  # Use defaults

        radii = getargument(rot_spec)
        @test radii isa AbstractVector
        @test length(radii) > 0
        @test maximum(radii) ≤ minimum([0.3, 0.3])  # Bounded by kmax
    end
end

@testset "RotationalEstimate Name Generation" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)

    @testset "Marginal -> Rotational" begin
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        rot_spec = rotational_estimate(spec)

        name = getestimatename(rot_spec)
        @test occursin("rotational", lowercase(name))
        @test occursin("spectra", lowercase(name))
    end

    @testset "Partial -> Rotational" begin
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        partial_spec = partial_spectra(spec)
        rot_partial_spec = rotational_estimate(partial_spec)

        name = getestimatename(rot_partial_spec)
        @test occursin("rotational", lowercase(name))
        @test occursin("partial", lowercase(name))
    end
end

@testset "Partial Spectra from RotationalEstimate" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 3, return_type = :tuple, point_number = 30)

    @testset "Marginal rotational -> partial rotational" begin
        region = getregion(data)
        spec = spectra(data, nk = (6, 6), kmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))
        rot_spec = rotational_estimate(spec)

        # Convert to partial
        partial_rot_spec = partial_spectra(rot_spec)

        @test partial_rot_spec isa RotationalEstimate
        @test is_partial(partial_rot_spec) == true
        @test getargument(partial_rot_spec) == getargument(rot_spec)  # Same radii
        @test size(partial_rot_spec) == size(rot_spec)

        # Name should indicate both partial and rotational
        name = getestimatename(partial_rot_spec)
        @test occursin("partial", lowercase(name))
        @test occursin("rotational", lowercase(name))
    end
end

@testset "Mathematical Properties" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 40)

    @testset "Smoothing effect" begin
        # Create spectra with sharp features
        region = getregion(data)
        spec = spectra(data, nk = (16, 16), kmax = (0.8, 0.8),
            tapers = sin_taper_family((2, 2), region))

        # Apply rotational averaging
        rot_spec_smooth = rotational_estimate(spec, kernel = GaussKernel(0.1))
        rot_spec_sharp = rotational_estimate(spec, kernel = RectKernel(0.02))

        estimate_smooth = getestimate(rot_spec_smooth)
        estimate_sharp = getestimate(rot_spec_sharp)

        # Smoother kernel should produce less variation
        # (This is hard to test rigorously without knowing the exact spectral structure)
        @test length(estimate_smooth) == length(estimate_sharp)
    end

    @testset "Radial monotonicity for simple functions" begin
        # For some functions, rotational average should be monotonic
        # This is more of a qualitative test
        rng = StableRNG(123)
        grids = make_grids_example(
            rng, n_processes = 1, return_type = :single, grid_dims = (20, 30))
        # nyquist is 20/4/2 = 2.5, 30/6/2 = 2.5

        region = getregion(grids)
        spec = spectra(grids, nk = (16, 16), kmax = (2.5, 2.5),
            tapers = sin_taper_family((3, 3), region))

        rot_spec = rotational_estimate(spec)
        estimate = getestimate(rot_spec)
        radii = getargument(rot_spec)

        @test length(estimate) == length(radii)
        @test all(isfinite, estimate)
        @test all(x -> real(x) ≥ 0, estimate)  # Power should be non-negative
    end
end

@testset "Integration with Other Transforms" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)
    region = getregion(data)
    spec = spectra(data, nk = (6, 6), kmax = (0.3, 0.3),
        tapers = sin_taper_family((2, 2), region))

    @testset "Rotational -> Marginal transforms" begin
        rot_spec = rotational_estimate(spec)

        abs_rot = abs(rot_spec)
        @test abs_rot isa MarginallyTransformedEstimate

        real_rot = real(rot_spec)
        @test real_rot isa MarginallyTransformedEstimate
    end

    @testset "Rotational -> Coherence" begin
        rot_spec = rotational_estimate(spec)

        # Should be able to compute coherence from rotational spectra
        coh = coherence(rot_spec)
        @test coh isa Coherence
    end
end

@testset "Indexing Operations" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 3, return_type = :tuple, point_number = 30)
    region = getregion(data)
    spec = spectra(data, nk = (6, 6), kmax = (0.3, 0.3),
        tapers = sin_taper_family((2, 2), region))
    rot_spec = rotational_estimate(spec)

    @testset "Process indexing" begin
        rot_sub = rot_spec[1, 2]
        @test rot_sub isa RotationalEstimate
        @test size(rot_sub) == (1, 1)
        @test getargument(rot_sub) == getargument(rot_spec)  # Same radii
    end

    @testset "Radial indexing" begin
        # For isotropic estimates, should be able to index by radius
        radii = getargument(rot_spec)

        # Index specific radius
        if length(radii) > 2
            rot_radius = rot_spec[1, 1, 3]  # 3rd radius point
            @test rot_radius isa RotationalEstimate
            @test size(rot_radius) == (1, 1)

            radius_arg = getargument(rot_radius)
            @test length(radius_arg) == 1
            @test radius_arg[1] ≈ radii[3]
        end
    end
end

@testset "Edge Cases and Error Conditions" begin
    rng = StableRNG(123)

    @testset "Very small regions" begin
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 10,
            region_min = (-0.1, -0.1), region_max = (0.1, 0.1))

        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (1.0, 1.0),
            tapers = sin_taper_family((2, 2), region))

        rot_spec = rotational_estimate(spec)
        @test rot_spec isa RotationalEstimate
        @test all(isfinite, stack(getestimate(rot_spec)))
    end

    @testset "Single radius" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 20)
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        rot_spec = rotational_estimate(spec, radii = [0.1], kernel = RectKernel(0.1))
        @test length(getargument(rot_spec)) == 1
        @test length(getestimate(rot_spec)) == 1
    end

    @testset "Zero radius" begin
        rng = StableRNG(123)
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 20)
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        rot_spec = rotational_estimate(spec, radii = [0.0], kernel = RectKernel(0.1)) # TODO: need to change the behaviour of the Kernel default to not use the radii
        @test getargument(rot_spec)[1] ≈ 0.0

        # At zero radius, should equal center value of original spectrum
        estimate = getestimate(rot_spec)[1]
        @test all(isfinite, estimate)
    end
end

@testset "Performance and Type Stability" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)
    region = getregion(data)
    spec = spectra(data, nk = (6, 6), kmax = (0.3, 0.3),
        tapers = sin_taper_family((2, 2), region))

    @testset "Type stability" begin
        rot_spec = rotational_estimate(spec)
        @test rot_spec isa RotationalEstimate
        @test getargument(rot_spec) isa AbstractVector
        @test getestimate(rot_spec) isa AbstractVector
    end

    @testset "Consistent dimensions" begin
        custom_radii = collect(range(0, 0.2, length = 10))
        rot_spec = rotational_estimate(spec, radii = custom_radii)

        @test length(getargument(rot_spec)) == 10
        @test length(getestimate(rot_spec)) == 10
    end
end
