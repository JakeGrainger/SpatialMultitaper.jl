@testset "apply_tapers" begin
    # TODO make a constant taper to check this
	tapers = interpolated_taper_family(
		[ones(5), ones(5)],
		CartesianGrid((-0.125,), (1.125,), dims = (5,)),
	)
	X = georef(
		(marks = ones(4),),
		PointSet([Point(0), Point(1), Point(0.5), Point(0.3)]),
	)
	@test SpatialMultitaper.apply_taper(domain(X), values(X)[1], tapers) ≈
		  fill(complex(1), 8) # 8 == ntapers * npoints
end