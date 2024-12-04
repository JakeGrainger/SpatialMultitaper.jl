@testset "sin taper" begin
	# check taper function is correct
	tapers = sin_taper_family((3, 3), Box(Point(0, 0), Point(10, 5)))
	@test tapers[6](0.4, 3.4) ≈
		  SpatialMultitaper.sin_taper(0.4, 3, 10) *
		  SpatialMultitaper.sin_taper(3.4, 2, 5)
	@test tapers[6](0.4, 3.4) ≈
		  sin(π * 3 * 0.4 / 10) * sin(π * 2 * 3.4 / 5) * 2 / sqrt(50)
	@test tapers[6](0.4, 3.4) ≈ tapers.tapers[6](0.4, 3.4)

	# TODO: check quality of interpolation methods

end

@testset "Taper" begin
	grid = CartesianGrid((0, 0), (10, 10), dims = (10, 10))
	disc_taper = SpatialMultitaper.discretetaper(ones(10, 10), grid)
	@test disc_taper isa SpatialMultitaper.DiscreteTaper

	int_taper = SpatialMultitaper.interpolate(disc_taper)
	@test int_taper isa SpatialMultitaper.InterpolatedTaper
	@test int_taper isa SpatialMultitaper.ContinuousTaper
	@test int_taper(20, 31.2) ≈ 0.0
	@test SpatialMultitaper.taper_ft(int_taper, (0.6, 3.4)) ≈ 0.0

	taper_family = interpolated_taper_family(
		[[zeros(5, 10); ones(5, 10)], [ones(5, 10); zeros(5, 10)]],
		grid,
	)
	@test taper_family[1](0.6, 3.4) ≈ 0.0
	@test taper_family[2](8, 3.4) ≈ 0.0
end
