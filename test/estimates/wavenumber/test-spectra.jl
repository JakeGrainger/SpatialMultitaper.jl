using SpatialMultitaper, Test, StableRNGs, StaticArrays
include("../../test_utilities/TestData.jl")
using .TestData

import SpatialMultitaper: Spectra, get_evaluation_points, get_estimates,
                          _dft_to_spectral_matrix!, SingleProcessTrait, MultipleVectorTrait,
                          MultipleTupleTrait, _compute_spectral_matrix,
                          get_estimation_information,
                          ProcessInformation, EstimationInformation, MarginalTrait,
                          MultipleVectorTrait, SingleProcessTrait, MultipleTupleTrait,
                          _compute_spectral_matrix!

@testset "Spectra Construction" begin
    rng = StableRNG(123)

    @testset "Basic construction" begin
        # Create test data
        wavenumber = (1:10, 1:10)
        power = rand(rng, ComplexF64, 2, 2, 10, 10)
        processinfo = ProcessInformation{2, MultipleVectorTrait}(
            [1, 2], [1, 2], ones(2, 2), ones(2, 2))
        estimationinfo = EstimationInformation(5)

        spec = Spectra{MarginalTrait}(wavenumber, power, processinfo, estimationinfo)
        @test get_evaluation_points(spec) == wavenumber
        @test get_estimates(spec) == power
        @test embeddim(spec) == 2
        @test size(spec) == (2, 2)
    end

    @testset "Single process case" begin
        wavenumber = (1:10,)
        power = rand(rng, Float64, 1, 1, 10)
        processinfo = ProcessInformation{1, MultipleVectorTrait}(
            [1], [1], ones(1, 1), ones(1, 1))
        estimationinfo = EstimationInformation(3)

        spec = Spectra{MarginalTrait}(wavenumber, power, processinfo, estimationinfo)
        @test size(spec) == (1, 1)
        @test embeddim(spec) == 1
    end
end

@testset "spectra function" begin
    rng = StableRNG(123)

    @testset "Points data" begin
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 50)
        nk = (8, 8)
        kmax = (0.5, 0.5)
        region = getregion(data)
        tapers = sin_taper_family((2, 2), region)

        spec = spectra(data, nk = nk, kmax = kmax, tapers = tapers)

        @test get_evaluation_points(spec) isa Tuple{AbstractVector, AbstractVector}
        @test length(get_evaluation_points(spec)[1]) == 8
        @test length(get_evaluation_points(spec)[2]) == 8
        @test embeddim(spec) == 2
        @test size(spec) == (2, 2)
        @test get_estimation_information(spec).ntapers == length(tapers)
    end

    @testset "Grid data" begin
        grids = make_grids_example(rng, n_processes = 1, return_type = :single,
            grid_dims = (10, 10), region_min = (0.0, 0.0), region_max = (10.0, 10.0))
        nk = (8, 8)
        kmax = (0.5, 0.5)  # Nyquist wavenumber for grid data
        region = getregion(grids)
        tapers = sin_taper_family((2, 2), region)

        spec = spectra(grids, nk = nk, kmax = kmax, tapers = tapers)
        @test size(spec) == ()
        @test embeddim(spec) == 2
    end

    @testset "Single process convenience" begin
        data = make_points_example(rng, n_processes = 1, return_type = :single)
        nk = (6, 6)
        kmax = (0.3, 0.3)
        region = getregion(data)
        tapers = sin_taper_family((2, 2), region)

        spec = spectra(data, nk = nk, kmax = kmax, tapers = tapers)
        @test size(spec) == ()
    end
end

@testset "DFT to Spectral Matrix" begin
    rng = StableRNG(123)

    @testset "MultipleVectorTrait" begin
        # P x M x n1 x n2 array
        J_n = rand(rng, ComplexF64, 3, 10, 8, 8)  # 3 processes, 10 tapers, 8x8 wavenumbers
        data = make_points_example(
            rng, n_processes = 3, return_type = :vector, point_number = 50)
        S_mat = zeros(ComplexF64, 3, 3, 8, 8)
        S_mat = _dft_to_spectral_matrix!(S_mat, J_n, MultipleVectorTrait())

        @test size(S_mat) == (3, 3, 8, 8)  # Should be spatial dimensions only
        @test eltype(S_mat) <: ComplexF64
    end

    @testset "SingleProcessTrait" begin
        # n1 x n2 x M array
        J_n = rand(rng, ComplexF64, 8, 8, 10)  # 8x8 wavenumbers, 10 tapers
        data = make_points_example(
            rng, n_processes = 1, return_type = :single, point_number = 50)
        S_mat = zeros(Float64, 8, 8)
        S_mat = _dft_to_spectral_matrix!(S_mat, J_n, SingleProcessTrait())

        @test size(S_mat) == (8, 8)
        @test eltype(S_mat) == Float64  # Should be real-valued for single process
    end

    @testset "MultipleTupleTrait" begin
        # Tuple of arrays n1 x n2 x M
        J_1 = rand(rng, ComplexF64, 8, 8, 10)
        J_2 = rand(rng, ComplexF64, 8, 8, 10)
        J_n = [SVector{2, ComplexF64}((J_1[i], J_2[i]))
               for i in CartesianIndices(size(J_1))]
        S_mat = Array{SMatrix{2, 2, ComplexF64, 4}}(undef, 8, 8)
        data = make_points_example(
            rng, n_processes = 2, return_type = :tuple, point_number = 50)

        S_mat = _dft_to_spectral_matrix!(S_mat, J_n, MultipleTupleTrait())
        @test size(S_mat) == (8, 8)
        @test eltype(S_mat) <: SMatrix{2, 2}

        # Single process case
        J_single = getindex.(J_n, Ref(SOneTo(1)))
        data_single = spatial_data((observations(data)[1],), getregion(data))
        S_single = Array{SMatrix{1, 1, ComplexF64, 1}}(undef, 8, 8)
        S_single = _dft_to_spectral_matrix!(S_single, J_single, MultipleTupleTrait())
        @test size(S_single) == (8, 8)
        @test eltype(S_single) <: SMatrix{1, 1}
    end
end

@testset "_compute_spectral_matrix function" begin
    rng = StableRNG(123)

    @testset "Vector input" begin
        x = rand(rng, ComplexF64, 5)
        S = _compute_spectral_matrix(x)
        @test S ≈ x * x'
        @test size(S) == (5, 5)
    end

    @testset "Matrix input" begin
        X = rand(rng, ComplexF64, 3, 10)  # 3 processes, 10 tapers
        S = zeros(ComplexF64, 3, 3)
        _compute_spectral_matrix!(S, X)
        expected = (X * X') / 10
        @test S ≈ expected
        @test size(S) == (3, 3)
    end
end

@testset "Spectra Indexing and Access" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 3, return_type = :tuple, point_number = 30)
    nk = (6, 6)
    kmax = (0.3, 0.3)
    region = getregion(data)
    tapers = sin_taper_family((2, 2), region)

    spec = spectra(data, nk = nk, kmax = kmax, tapers = tapers)

    @testset "Process indexing" begin
        spec_sub = spec[1, 2]
        @test size(spec_sub) == ()  # Single process pair
        @test get_evaluation_points(spec_sub) == get_evaluation_points(spec)  # Same wavenumbers
    end

    @testset "Wavenumber indexing" begin
        spec_wavenumber = spec[1, 1, 3, 4]  # specific wavenumber bin
        @test size(spec_wavenumber) == ()
        wavenumber_arg = get_evaluation_points(spec_wavenumber)
        @test length(wavenumber_arg[1]) == 1  # Single wavenumber
        @test length(wavenumber_arg[2]) == 1
    end
end

@testset "Edge Cases and Error Handling" begin
    rng = StableRNG(123)

    @testset "Empty or minimal data" begin
        # Test with very small datasets
        data = spatial_data(PointSet([Point(0.0, 0.0)]), Box(Point(-1, -1), Point(1, 1)))
        nk = (2, 2)
        kmax = (0.1, 0.1)
        region = getregion(data)
        tapers = sin_taper_family((1, 1), region)

        # Should not crash but might have limited meaningful results
        spec = spectra(data, nk = nk, kmax = kmax, tapers = tapers)
        @test size(spec) == ()
    end
end

@testset "automatic tapers" begin
    rng = StableRNG(123)
    data = make_points_example(rng, n_processes = 1, return_type = :single)
    nk = (6, 6)
    kmax = (0.3, 0.3)
    region = getregion(data)
    tapers = sin_taper_family((3, 3), region)
    spec = spectra(data, nk = nk, kmax = kmax, nw = 2) # so number of tapers = 2*2-1=3
    spec2 = spectra(data, nk = nk, kmax = kmax, tapers = tapers)
    @test spec == spec2
end

# TODO: need to add this
# @testset "NaN handling" begin
