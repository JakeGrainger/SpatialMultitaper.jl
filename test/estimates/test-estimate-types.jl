using SpatialMultitaper, Test
include("../test_utilities/TestUtils.jl")
using .TestUtils

import SpatialMultitaper: EstimateTrait, MarginalTrait, PartialTrait,
                          AbstractEstimate, AnisotropicEstimate, IsotropicEstimate,
                          ProcessInformation, EstimationInformation,
                          is_partial, checkprocessinformation, checkinputs,
                          SingleProcessTrait, MultipleVectorTrait, MultipleTupleTrait

@testset "Traits and Type System" begin
    @testset "EstimateTrait hierarchy" begin
        @test MarginalTrait <: EstimateTrait
        @test PartialTrait <: EstimateTrait
    end

    @testset "Abstract type hierarchy" begin
        @test AnisotropicEstimate{MarginalTrait, 2, 1, 1} <:
              AbstractEstimate{MarginalTrait, 2, 1, 1, 2}
        @test IsotropicEstimate{MarginalTrait, 2, 1, 1} <:
              AbstractEstimate{MarginalTrait, 2, 1, 1, 1}
    end
end

@testset "ProcessInformation" begin
    @testset "Construction" begin
        # Test 1D case
        processinfo_1d = ProcessInformation{1, SingleProcessTrait}(1, 1, 2, 2)
        @test processinfo_1d.process_indices_1 == 1
        @test processinfo_1d.process_indices_2 == 1
        @test processinfo_1d.mean_product == 2
        @test embeddim(processinfo_1d) == 1

        # Test 2D case
        processinfo_2d = ProcessInformation{2, MultipleVectorTrait}(
            [1, 2], [1, 2], ones(2, 2), ones(2, 2))
        @test length(processinfo_2d.process_indices_1) == 2
        @test embeddim(processinfo_2d) == 2
    end

    @testset "checkprocessinformation" begin
        # Valid cases
        processinfo = ProcessInformation{2, MultipleVectorTrait}(
            [1, 2], [1, 2], ones(2, 2), ones(2, 2))
        @test checkprocessinformation(processinfo, randn(2, 2, 10, 10)) == (2, 2)

        # Invalid cases
        @test_throws ArgumentError checkprocessinformation(
            processinfo, randn(3, 2, 10, 10))  # Wrong P
        @test_throws ArgumentError checkprocessinformation(
            processinfo, randn(2, 3, 10, 10))  # Wrong Q

        @test_throws ArgumentError ProcessInformation{2, MultipleVectorTrait}(
            [1, 2], [1, 2], ones(2, 3), ones(2, 2))  # Wrong mean_product size
    end
end

@testset "EstimationInformation" begin
    @test EstimationInformation(10).ntapers == 10
    @test EstimationInformation(nothing).ntapers === nothing
end

@testset "Bounds Checking" begin
    # Create mock estimate type for testing
    struct MockEstimate{E, D, P, Q, N} <: AbstractEstimate{E, D, P, Q, N}
        argument::NTuple{N}
        estimate::AbstractArray
        processinformation::ProcessInformation
        estimationinformation::EstimationInformation
        function MockEstimate{E}(
                argument::NTuple{N}, estimate, processinfo::ProcessInformation{D},
                estimationinfo) where {E, D, N}
            P, Q = checkinputs(argument, estimate, processinfo)
            new{E, D, P, Q, N}(argument, estimate, processinfo, estimationinfo)
        end
    end
    SpatialMultitaper.getargument(est::MockEstimate) = est.argument
    SpatialMultitaper.getestimate(est::MockEstimate) = est.estimate

    # Test with matrix estimate
    mock_matrix = MockEstimate{MarginalTrait}(
        (1:10, 1:10),
        rand(2, 2, 10, 10),
        ProcessInformation{2, MultipleVectorTrait}([1, 2], [1, 2], ones(2, 2), ones(2, 2)),
        EstimationInformation(5)
    )

    @testset "Process bounds checking" begin
        @test checkbounds(mock_matrix, 1, 1) === nothing
        @test checkbounds(mock_matrix, 2, 2) === nothing
        @test_throws BoundsError checkbounds(mock_matrix, 3, 1)
        @test_throws BoundsError checkbounds(mock_matrix, 1, 3)
    end

    @testset "Index bounds checking" begin
        @test checkbounds(mock_matrix, 1, 1, 5, 5) === nothing
        @test_throws BoundsError checkbounds(mock_matrix, 1, 1, 11, 5)
        @test_throws BoundsError checkbounds(mock_matrix, 1, 1, 5, 11)
    end
end

@testset "Input Validation" begin
    @testset "Array + ProcessInformation validation" begin
        # Valid 1D case
        processinfo = ProcessInformation{1, MultipleVectorTrait}(
            [1], [1], ones(1, 1), ones(1, 1))
        @test checkinputs((1:10,), rand(1, 1, 10), processinfo) == (1, 1)

        # Valid 2D case with matrices
        pi_2d = ProcessInformation{2, MultipleVectorTrait}(
            [1, 2], [1, 2], ones(2, 2), ones(2, 2))
        @test checkinputs((1:10, 1:10), rand(2, 2, 10, 10), pi_2d) == (2, 2)

        # Invalid dimension mismatch
        @test_throws ArgumentError checkinputs((1:10,), rand(2, 2, 10, 10), pi_2d)  # Wrong dimensions
        @test_throws ArgumentError checkinputs((1:10, 1:5), rand(2, 2, 10, 10), pi_2d)  # Size mismatch
    end
end

@testset "Utility Functions" begin
    @testset "is_partial" begin
        # Would test with actual estimate instances
        # @test is_partial(marginal_estimate) == false
        # @test is_partial(partial_estimate) == true
    end

    @testset "Type information" begin
        mock = MockEstimate{MarginalTrait}(
            (1:10, 1:10), rand(2, 2, 10, 10),
            ProcessInformation{2, MultipleVectorTrait}(
                [1, 2], [1, 2], ones(2, 2), ones(2, 2)),
            EstimationInformation(5)
        )

        @test embeddim(mock) == 2
        @test size(mock) == (2, 2)
    end
end
