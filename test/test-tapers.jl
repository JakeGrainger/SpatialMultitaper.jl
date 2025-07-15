@testset "sin taper" begin
    # check taper function is correct
    tapers = sin_taper_family((3, 3), Box(Point(0, 0), Point(10, 5)))
    @test tapers[6](0.4, 3.4) ≈ Spmt.sin_taper(0.4, 3, 10) * Spmt.sin_taper(3.4, 2, 5)
    @test tapers[6](0.4, 3.4) ≈ sin(π * 3 * 0.4 / 10) * sin(π * 2 * 3.4 / 5) * 2 / sqrt(50)
    @test tapers[6](0.4, 3.4) ≈ tapers.tapers[6](0.4, 3.4)

    tapers_alt = sin_taper_family((3, 3), Box(Point(2.3, 0), Point(12.3, 5)))
    @test tapers_alt[3](2.7, 1.4) ≈ tapers[3](0.4, 1.4)
end

@testset "interpolation" begin
    # TODO: check quality of interpolation methods
    x = ones(100)
    g = CartesianGrid((0,), (100,), dims = (100,))
    @test Spmt.L2_inner_product_interpolated(x, x, g) ≈ 99 + 2 / 3
end
@testset "Taper" begin
    grid = CartesianGrid((0, 0), (10, 10), dims = (10, 10))
    disc_taper = Spmt.discretetaper(ones(10, 10), grid)
    @test disc_taper isa Spmt.DiscreteTaper

    taper_family = interpolated_taper_family(
        [[zeros(5, 10); ones(5, 10)], [ones(5, 10); zeros(5, 10)]],
        grid,
    )
    @test taper_family[1](0.6, 3.4) ≈ 0.0
    @test taper_family[2](8, 3.4) ≈ 0.0

    int_taper = taper_family[1]
    @test int_taper isa Spmt.InterpolatedTaper
    @test int_taper isa Spmt.ContinuousTaper
    @test int_taper(20, 31.2) ≈ 0.0
    @test Spmt.taper_ft(int_taper, (0.6, 3.4)) ≈ 0.0
end
