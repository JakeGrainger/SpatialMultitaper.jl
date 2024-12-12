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

		dft = Spmt.tapered_dft(
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