region = Box(Point(0, 0), Point(3, 3))
pattern = georef((marks = [1, 2.4],), PointSet([Point(0, 0), Point(1, 1)]))
pattern2 = PointSet([Point(0, 0.3), Point(0.2, 0.2), Point(0.1, 1)])
griddata = georef((rf = rand(6 * 6),), CartesianGrid((0, 0), (3, 3), dims = (6, 6)))
tapers = sin_taper_family((3, 3), region)
nfreq = (10, 10)
fmax = (2, 2)
mt_est = multitaper_estimate(
	(pattern, pattern2, griddata),
	region,
	nfreq = nfreq,
	fmax = fmax,
	tapers = tapers,
)

x = complex_coherence(mt_est)
@test x isa NamedTuple
@test x.freq == mt_est.freq
@test x.complex_coherence isa Array{<:Number, 4}
@test size(x.complex_coherence) == size(mt_est.power)
@test x.complex_coherence_jackknifed === nothing

x = magnitude_coherence(mt_est)
@test x isa NamedTuple
@test x.freq == mt_est.freq
@test x.magnitude_coherence isa Array{<:Real, 4}
@test size(x.magnitude_coherence) == size(mt_est.power)
@test x.magnitude_coherence_jackknifed === nothing

x = group_delay(mt_est)
@test x isa NamedTuple
@test x.freq == mt_est.freq
@test x.group_delay isa Array{<:Real, 4}
@test size(x.group_delay) == size(mt_est.power)
@test x.group_delay_jackknifed === nothing

x = partial_complex_coherence(mt_est)
@test x isa NamedTuple
@test x.freq == mt_est.freq
@test x.partial_complex_coherence isa Array{<:Number, 4}
@test size(x.partial_complex_coherence) == size(mt_est.power)
@test x.partial_complex_coherence_jackknifed === nothing

x = partial_magnitude_coherence(mt_est)
@test x isa NamedTuple
@test x.freq == mt_est.freq
@test x.partial_magnitude_coherence isa Array{<:Real, 4}
@test size(x.partial_magnitude_coherence) == size(mt_est.power)
@test x.partial_magnitude_coherence_jackknifed === nothing

x = partial_group_delay(mt_est)
@test x isa NamedTuple
@test x.freq == mt_est.freq
@test x.partial_group_delay isa Array{<:Real, 4}
@test size(x.partial_group_delay) == size(mt_est.power)
@test x.partial_group_delay_jackknifed === nothing

x = partial_spectra(mt_est)
x_alt = partial_spectra(mt_est, 1, 1, 2:3, 2:3)
@test x.partial_spectra[1, 1, :, :] ≈ x_alt.partial_spectra[1, 1, :, :]

function long_partial_spectra(x)
	C = inv(x)
	par_coh = -complex_coherence(C)
	Sa = inv.(SpatialMultitaper.diag(C))
	return [
		i == j ? Sa[i] : par_coh[i, j] / (1 - abs2(par_coh[i, j])) * sqrt(Sa[i] * Sa[j]) for
		i in axes(C, 1), j in axes(C, 2)
	]
end
@test long_partial_spectra(mt_est.power[:,:,1,1]) ≈ x.partial_spectra[:,:,1,1]