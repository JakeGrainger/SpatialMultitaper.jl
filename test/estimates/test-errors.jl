using SpatialMultitaper, Test
include("../test_utilities/TestUtils.jl")
using .TestUtils

import SpatialMultitaper: partial_from_marginal_error, KFunction, CFunction

@testset "Error Messages" begin
    @testset "partial_from_marginal_error" begin
        # Test the error message generation
        error_msg = partial_from_marginal_error(KFunction, CFunction)

        @test error_msg isa ArgumentError
        @test occursin("Cannot compute partial", string(error_msg))
        @test occursin("K function", string(error_msg))  # From KFunction base name
        @test occursin("C function", string(error_msg))  # From CFunction base name
        @test occursin("k_function", string(error_msg))  # Suggests alternative
    end

    @testset "Error message contains helpful suggestions" begin
        error_msg = partial_from_marginal_error(KFunction, CFunction)
        error_str = string(error_msg)

        # Should suggest alternatives
        @test occursin("original spectrum", error_str) ||
              occursin("original data", error_str)
        @test occursin("k_function", error_str)  # Function name suggestion
    end
end

@testset "Integration with Actual Error Scenarios" begin
    # Test that the error is actually thrown in appropriate contexts
    # This would depend on how the error checking is implemented in the actual functions

    @testset "Placeholder for actual error scenarios" begin
        # When the actual spatial functions (K, L, C functions) are implemented,
        # this is where you'd test that trying to compute partial versions
        # from marginal spatial functions throws the appropriate error

        @test_skip "Error throwing integration tests - implement when spatial functions are complete"
    end
end
