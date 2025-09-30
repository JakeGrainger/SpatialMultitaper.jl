using SpatialMultitaper, Test, StableRNGs
include("../test_utilities/TestData.jl")
using .TestData

@testset "Estimate Printing" begin
    rng = StableRNG(123)
    data = make_points_example(
        rng, n_processes = 2, return_type = :tuple, point_number = 30)

    @testset "Spectra printing" begin
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))

        # Test compact show
        io = IOBuffer()
        show(io, spec)
        output = String(take!(io))

        @test occursin("spectra", lowercase(output))
        @test occursin("2D", output)  # Dimension
        @test occursin("↔", output)   # Process separator

        # Test detailed show
        io = IOBuffer()
        show(io, MIME"text/plain"(), spec)
        detailed_output = String(take!(io))

        @test occursin("dimensional process", detailed_output)
        @test occursin("between processes", detailed_output)
        @test occursin("evaluated at", detailed_output)
        @test occursin("with values of type", detailed_output)
    end

    @testset "Coherence printing" begin
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        coh = coherence(spec)

        io = IOBuffer()
        show(io, coh)
        output = String(take!(io))

        @test occursin("coherence", lowercase(output))
        @test occursin("2D", output)
    end

    @testset "Partial estimates printing" begin
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        partial_spec = partial_spectra(spec)

        io = IOBuffer()
        show(io, partial_spec)
        output = String(take!(io))

        @test occursin("partial", lowercase(output))
        @test occursin("spectra", lowercase(output))
    end

    @testset "Transformed estimates printing" begin
        region = getregion(data)
        spec = spectra(data, nk = (4, 4), kmax = (0.2, 0.2),
            tapers = sin_taper_family((2, 2), region))
        abs_spec = abs(spec)

        io = IOBuffer()
        show(io, abs_spec)
        output = String(take!(io))

        @test occursin("abs", lowercase(output))
        @test occursin("spectra", lowercase(output))
    end

    @testset "Rotational estimates printing" begin
        region = getregion(data)
        spec = spectra(data, nk = (6, 6), kmax = (0.3, 0.3),
            tapers = sin_taper_family((2, 2), region))
        rot_spec = rotational_estimate(spec)

        io = IOBuffer()
        show(io, rot_spec)
        output = String(take!(io))

        @test occursin("rotational", lowercase(output))
        @test occursin("spectra", lowercase(output))
    end
end

@testset "Argument Printing" begin
    @testset "_printargument with tuples" begin
        import SpatialMultitaper: _printargument

        # Test tuple printing
        arg_tuple = (1:5, 2:6)
        output = _printargument(arg_tuple)
        @test occursin(",", output)
        @test occursin(":", output)  # Range notation
    end

    @testset "_printargument with single values" begin
        import SpatialMultitaper: _printargument

        single_arg = 1:10
        output = _printargument(single_arg)
        @test isa(output, String)
        @test length(output) > 0
    end
end

@testset "Process Names Display" begin
    rng = StableRNG(123)

    @testset "Multiple processes" begin
        data = make_points_example(
            rng, n_processes = 3, return_type = :tuple, point_number = 20)
        region = getregion(data)
        spec = spectra(data, nk = (3, 3), kmax = (0.15, 0.15),
            tapers = sin_taper_family((2, 2), region))

        io = IOBuffer()
        show(io, MIME"text/plain"(), spec)
        output = String(take!(io))

        # Should show process indices
        @test occursin("[1, 2, 3]", output) || occursin("1:3", output)
    end

    @testset "Single process" begin
        data = make_points_example(
            rng, n_processes = 1, return_type = :single, point_number = 20)
        region = getregion(data)
        spec = spectra(data, nk = (3, 3), kmax = (0.15, 0.15),
            tapers = sin_taper_family((2, 2), region))

        io = IOBuffer()
        show(io, MIME"text/plain"(), spec)
        output = String(take!(io))

        @test occursin("[1]", output) || occursin("1", output)
    end
end

@testset "Different Dimensions Display" begin
    rng = StableRNG(123)

    @testset "1D case" begin
        # Create 1D data
        data = make_points_example(rng, n_processes = 2, return_type = :tuple,
            point_number = 20, dim = 1, region_min = (0,), region_max = (1,))

        region = getregion(data)
        spec = spectra(data, nk = (10,), kmax = (0.5,),
            tapers = sin_taper_family((3,), region))

        io = IOBuffer()
        show(io, spec)
        output = String(take!(io))

        @test occursin("1D", output)
    end
end

@testset "Edge Cases in Printing" begin
    @testset "Very long process lists" begin
        # Test with many processes to see if output is reasonable
        rng = StableRNG(123)
        data = make_grids_example(rng, n_processes = 10, return_type = :vector,
            grid_dims = (5, 5), region_min = (0.0, 0.0), region_max = (1.0, 1.0))

        region = getregion(data)
        spec = spectra(data, nk = (11, 11), kmax = (2.5, 2.5),
            tapers = sin_taper_family((2, 2), region))

        io = IOBuffer()
        show(io, spec)
        output = String(take!(io))

        # Should not be unreasonably long
        @test length(output) < 1000  # Reasonable upper bound
        @test occursin("10", output) # Should mention the number of processes
    end
end
