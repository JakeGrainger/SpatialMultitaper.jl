abstract type SpatialShift end
struct UniformShift{D, T} <: SpatialShift
	min::NTuple{D, T}
	max::NTuple{D, T}
end
Base.rand(u::UniformShift) = rand.(Uniform.(u.min, u.max))

struct UniformBallShift{D,T<:Number} <: SpatialShift
	max_shift::T
	UniformBallShift(max_shift::T, dim::Val{D}) where {T, D} = new{D, T}(max_shift)
end
function Base.rand(u::UniformShift{2})
	radius = random(Uniform(0, u.max_shift))
	angle = random(Uniform(0, 2π))
	return (radius * cos(angle), radius * sin(angle))
end
function Base.rand(u::UniformShift{3})
	radius = random(Uniform(0, u.max_shift))
	theta = random(Uniform(0, π))
	phi = random(Uniform(0, 2π))
	return (
		radius * sin(theta) * cos(phi),
		radius * sin(theta) * sin(phi),
		radius * cos(theta),
	)
end
Base.rand(u::UniformBallShift{D,T}) where {D,T} = error("Unsupported dimension $D for UniformBallShift")


##
abstract type ShiftMethod end
struct NoShift <: ShiftMethod end
marginal_shift(pp::PointSet, ::NoShift) = pp

struct ToroidalShift{R <: Box, S <: Union{<:SpatialShift, <:NTuple}} <: ShiftMethod
	region::R
	shift::S
end
function ToroidalShift(box::Box)
	centered_box = inverse(Translate(to(centroid(box))...))(box)
	ToroidalShift(
		box,
		UniformShift(unitless_coords(centered_box.min), unitless_coords(centered_box.max)),
	)
end
Base.rand(shift::ToroidalShift) = ToroidalShift(shift.region, rand(shift.shift))

function toroidal_shift(pp::PointSet, region::Box, shift)
	return PointSet([toroidal_shift(p, region, shift) for p in pp])
end

function toroidal_shift(p::Point, region::Box, shift::Tuple)
	sides = SpatialMultitaper.box2sides(region)
	return Point(toroidal_shift.(SpatialMultitaper.unitless_coords(p), sides, shift))
end

function toroidal_shift(x, side, v)
	a = side[1]
	b = side[2]
	return a + mod(x - a + v, b - a)
end
marginal_shift(pp::PointSet, shift_method::ToroidalShift) =
	toroidal_shift(pp, shift_method.region, shift_method.shift)

struct MinusShift{R,G,S}
	region::R
	inset_region::G
	shift::S
end
function MinusShift(region, maxshift)
	inset_region = make_inset_region(region, maxshift)
	shift = UniformBallShift(maxshift, Val{embeddim(region)}())
	return MinusShift(region, inset_region, shift)
end

function make_inset_region(region, maxshift)
	bbox = boundingbox(region)
	bsize = map(x->x[2]-x[1], box2sides(bbox))
	transform = Stretch((1 .- 2maxshift./bsize)...)
	return transform(region)
end

##
function shift_resample(
	data::NTuple{P, S},
	region,
	groups,
	statistic,
	shift_method::ShiftMethod,
) where {P, S}
	@assert sort(reduce(vcat, groups)) == 1:P "groups of shifts should partition the space"
	group_shifts = Dict(group => rand(shift_method) for group in groups)
	shifted_processes =
		ntuple(p -> marginal_shift(data[p], group_shifts[findgroup(p, groups)]), Val{P}())
	statistic(shifted_processes, region)
end

findgroup(p, groups) = groups[findfirst(g -> p ∈ g, groups)]

##
function partial_shift_resample(
	data::NTuple{P, S},
	region,
	statistic,
	shift_method::ShiftMethod,
) where {P, S}
	groups = [1:P, P+1:2P]
	augmented_data = (data..., deepcopy(data)...)
	return shift_resample(augmented_data, region, groups, statistic, shift_method)
end

function partial_K_resample(
	data::NTuple,
	region;
	radii,
	shift_method::ShiftMethod,
	tapers,
	nfreq,
	fmax,
)
	p = length(data)
	firsthalf = 1:p
	secondhalf = p+1:2p
	indices = [
		(x, y, view(firsthalf, Not(SVector(i, j))), view(secondhalf, Not(SVector(i, j)))) for (i, x) in enumerate(firsthalf), (j, y) in enumerate(secondhalf) if i <= j
	]
	function wrapped_partial_K(_data, _region)
		partial_K(
			_data,
			_region;
			radii = radii,
			tapers = tapers,
			nfreq = nfreq,
			fmax = fmax,
			indices = indices,
		)
	end
	resampled = partial_shift_resample(data, region, wrapped_partial_K, shift_method)
	return (
		radii = resampled.radii,
		partial_K = Dict((key[1], key[2] - p) => val for (key, val) in resampled.partial_K),
	)
end

## freq domain null resampling
function multitaper_estimate_resampled(
    data,
    region;
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
	nresamples::Int
)
    data, dim = check_spatial_data(data)
    mean_method = check_mean_method(mean_method, data)
    freq = make_freq(nfreq, fmax, dim)
    J_n = tapered_dft(data, tapers, nfreq, fmax, region, mean_method)
	power = dft2spectralmatrix(J_n)
	resampled = [SpectralEstimate(freq, dft2spectralmatrix(null_resample(J_n))) for _ in 1:nresamples]
	observed = SpectralEstimate(freq, power)
    return  (observed = observed, resampled = resampled)
end

"""
	null_resample(J_n::Array)

Resamples the DFTs across M independently in P, assuming that the DFTs are stored as one large array which is P x M x n_1 x ... x n_D.
"""
function null_resample(J_n::Array{T,N}) where {T,N}
	error("under development!!")
	M = size(J_n, 2)
	P = size(J_n, 1)
	permutedims(map(j -> selectdim(selectdim(J_n, 1, j), 2, rand(1:M, M)), 1:P), (N, 1:N-1...))
end

"""
	null_resample(J_n::NTuple{P, Array{T, N}}) where {P, T, N}

Resamples the DFTs across M independently in P, assuming that the DFTs are stored as a tuple of P arrays of size n_1 x ... x n_D x M.
"""
function null_resample(J_n::NTuple{P,Array{T,N}}) where {P,T,N}
	M = size(first(J_n))[end]
	# ntuple(j -> selectdim(J_n[j], N, rand(1:M, M)), Val{P}())
	new_J = ntuple(p -> Array{T,N}(undef, size(J_n[p])), Val{P}())
	for p in 1:P
		for i in CartesianIndices(J_n[p])
			new_J[p][i] = J_n[p][i.I[1:end-1]..., rand(1:M)]
		end
	end
	return new_J
end