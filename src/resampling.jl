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

struct ToroidalShift{R<:Box,S<:Union{<:SpatialShift, <:NTuple}} <: ShiftMethod
    region::R
    shift::S
end
function ToroidalShift(box::Box)
    centered_box = inverse(Translate(to(centroid(box))...))(box)
    ToroidalShift(box, UniformShift(centered_box.min, centered_box.max))
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
marginal_shift(pp::PointSet, shift_method::ToroidalShift) = toroidal_shift(pp, shift_method.region, shift_method.shift)

##
function shift_resample(data::NTuple{P,S}, groups, region, statistic, shift_method::ShiftMethod) where {P,S}
    @assert sort(reduce(vcat,groups)) == 1:P "groups of shifts should partition the space"
    group_shifts = Dict(group => rand(shift) for group in groups)
    shifted_processes = (p->marginal_shift(data[p], group_shifts[findgroup[p, groups]]), Val{P}())
    statistic(shifted_processes, region)
end

findgroup(p, groups) = findfirst(g->p ∈ g, groups)