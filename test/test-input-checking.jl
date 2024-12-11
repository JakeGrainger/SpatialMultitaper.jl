@testset "check_spatial_data" begin
	p1 = georef((rf = ones(2),), PointSet([Point(0), Point(1)]))
	p2 = georef((rf = ones(2),), PointSet([Point(0, 0), Point(1, 1)]))
	p3 = PointSet([Point(0, 0), Point(1, 1)])
	@test ((p2, p3), 2) == SpatialMultitaper.check_spatial_data((p2, p3))
	@test_throws AssertionError SpatialMultitaper.check_spatial_data((p1, p2))
	@test_throws AssertionError SpatialMultitaper.check_spatial_data((p1, p3))

	pattern = georef((marks = [1],), PointSet([Point(0), Point(1)]))
	griddata = georef((rf = rand(6 * 6),), CartesianGrid((0, 0), (3, 3), dims = (6, 6)))
	multi_grid = georef(
		(rf1 = griddata.rf, rf2 = rand(6 * 6)),
		CartesianGrid((0, 0), (3, 3), dims = (6, 6)),
	)

	@test_throws AssertionError SpatialMultitaper.check_spatial_data((
		pattern,
		griddata,
	))
	@test_throws AssertionError SpatialMultitaper.check_spatial_data((
		domain(pattern),
		griddata,
	))
	@test_warn "more than one random field provided to a geotable, currently we only process the first of these!" SpatialMultitaper.check_spatial_data((
		multi_grid,
	))

end
@testset "jk_weight" begin
	x = rand(30, 4)
	@test SpatialMultitaper.spectral_matrix(x[[1:20; 22:30], :]) ==
		  SpatialMultitaper.spectral_matrix(x, SpatialMultitaper.make_jk_weight(30, 21))
end
@testset "apply_tapers" begin
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
@testset "pipeline" begin
	@testset "1d" begin
		region = Box(Point(0), Point(3))
		pattern = georef(
			(a = ones(4),),
			PointSet([Point(0), Point(1), Point(0.5), Point(0.3)]),
		)
		griddata = georef((rf = rand(6),), CartesianGrid((0,), (3,), dims = (6,)))
		tapers = interpolated_taper_family(
			[ones(3), ones(3), ones(3)],
			CartesianGrid((0,), (3,), dims = (3,)),
		)
		nfreq = (10,)
		fmax = (2,)
		mt_est = multitaper_estimate(
			(pattern, griddata),
			region,
			nfreq = nfreq,
			fmax = fmax,
			tapers = tapers,
		)
		@test size(mt_est.power) == (2, 2, 10)
		mt_est_jk = multitaper_estimate(
			(pattern, griddata),
			region,
			nfreq = nfreq,
			fmax = fmax,
			tapers = tapers,
			jackknife = true,
		)
		@test size(mt_est_jk.power_jackknifed) == (3,)
		all(size.(mt_est_jk.power_jackknifed) .== (2, 2, 10))
		@test true

		mt_est_1 = multitaper_estimate(
			pattern,
			region,
			nfreq = nfreq,
			fmax = fmax,
			tapers = tapers,
		)
		@test mt_est_1.power == multitaper_estimate(
			(pattern,),
			region,
			nfreq = nfreq,
			fmax = fmax,
			tapers = tapers,
		).power[
			1,
			1,
			:,
		]
		@test size(mt_est_1.power) == (10,)
	end
	@testset "2d" begin
		region = Box(Point(0, 0), Point(3, 3))
		pattern = georef((marks = [1, 2.4],), PointSet([Point(0, 0), Point(1, 1)]))
		griddata =
			georef((rf = rand(6 * 6),), CartesianGrid((0, 0), (3, 3), dims = (6, 6)))
		tapers = interpolated_taper_family(
			[ones(3, 3), ones(3, 3), ones(3, 3)],
			CartesianGrid((0, 0), (3, 3), dims = (3, 3)),
		)
		nfreq = (10, 10)
		fmax = (2, 2)
		mt_est = multitaper_estimate(
			(pattern, griddata),
			region,
			nfreq = nfreq,
			fmax = fmax,
			tapers = tapers,
		)
		@test size(mt_est.power) == (2, 2, 10, 10)

		dft = SpatialMultitaper.tapered_dft(
			(pattern, griddata),
			tapers,
			nfreq,
			fmax,
			region,
			(DefaultMean(), DefaultMean()),
		)
		@test mt_est.power[1, 2, 5, 2] ≈
			  sum(dft[:, 1, 5, 2] .* conj.(dft[:, 2, 5, 2])) / 3
	end
end
@testset "NaN handling" begin

	griddata = georef((rf = rand(6 * 6),), CartesianGrid((0, 0), (3, 3), dims = (6, 6)))
	multi_grid = georef(
		(rf1 = griddata.rf, rf2 = rand(6 * 6)),
		CartesianGrid((0, 0), (3, 3), dims = (6, 6)),
	)
	nan_grid =
		georef((rf = fill(NaN, 6 * 6),), CartesianGrid((0, 0), (3, 3), dims = (6, 6)))

	region = Box(Point(0, 0), Point(3, 3))
	tapers = interpolated_taper_family(
		[ones(3, 3), ones(3, 3), ones(3, 3)],
		CartesianGrid((0, 0), (3, 3), dims = (3, 3)),
	)
	nfreq = (10, 10)
	fmax = (2, 2)

	@test multitaper_estimate(
		griddata,
		region,
		nfreq = nfreq,
		fmax = fmax,
		tapers = tapers,
	).power ≈
		  multitaper_estimate(
		multi_grid,
		region,
		nfreq = nfreq,
		fmax = fmax,
		tapers = tapers,
	).power
	@test multitaper_estimate(
		nan_grid,
		region,
		nfreq = nfreq,
		fmax = fmax,
		tapers = tapers,
	).power == zeros(10, 10)
end
