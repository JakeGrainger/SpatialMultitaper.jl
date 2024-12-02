abstract type SpatialShift end
struct UniformShift{D,T} <: SpatialShift
    min::NTuple{D,T}
    max::NTuple{D,T}
end
Base.rand(u::UniformShift) = rand.(Uniform.(u.min, u.max))

##
abstract type ShiftMethod end
struct NoShift <: ShiftMethod end
marginal_shift(pp::PointSet, ::NoShift) = pp

struct ToroidalShift{R<:Box,S<:Union{<:SpatialShift,<:NTuple}} <: ShiftMethod
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

##
function shift_resample(
    data::NTuple{P,S},
    groups,
    region,
    statistic,
    shift_method::ShiftMethod,
) where {P,S}
    @assert sort(reduce(vcat, groups)) == 1:P "groups of shifts should partition the space"
    group_shifts = Dict(group => rand(shift_method) for group in groups)
    shifted_processes =
        ntuple(p -> marginal_shift(data[p], group_shifts[findgroup(p, groups)]), Val{P}())
    statistic(shifted_processes, region)
end

findgroup(p, groups) = groups[findfirst(g -> p ∈ g, groups)]

##
function partial_shift_resample(
    data::NTuple{P,S},
    region,
    statistic,
    shift_method::ShiftMethod,
) where {P,S}
    groups = [1:P, P+1:2P]
    augmented_data = (data..., deepcopy(data)...)
    return shift_resample(augmented_data, groups, region, statistic, shift_method)
end

function partial_K_resample(
    data::NTuple,
    radii,
    tapers,
    shift_method::ShiftMethod;
    region,
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
        partial_K(_data, radii, tapers; region = _region, nfreq = nfreq, fmax = fmax, indices = indices)
    end
    resampled = partial_shift_resample(data, region, wrapped_partial_K, shift_method)
    return (radii = resampled.radii, partial_K = Dict((key[1], key[2] - p) => val for (key, val) in resampled.partial_K))
end