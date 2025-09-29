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
