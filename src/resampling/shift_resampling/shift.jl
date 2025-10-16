abstract type SpatialShift end
struct UniformShift{D, T} <: SpatialShift
    min::NTuple{D, T}
    max::NTuple{D, T}
end
Base.rand(rng::AbstractRNG, u::UniformShift) = rand.(rng, Uniform.(u.min, u.max))

struct UniformBallShift{D, T <: Number} <: SpatialShift
    max_shift::T
    UniformBallShift(max_shift::T, dim::Val{D}) where {T, D} = new{D, T}(max_shift)
end
function Base.rand(rng::AbstractRNG, u::UniformBallShift{1})
    return rand(rng, Uniform(-u.max_shift, u.max_shift))
end
function Base.rand(rng::AbstractRNG, u::UniformBallShift{2})
    radius = rand(rng, Uniform(0, u.max_shift))
    angle = rand(rng, Uniform(0, 2π))
    return (radius * cos(angle), radius * sin(angle))
end
function Base.rand(rng::AbstractRNG, u::UniformBallShift{3})
    radius = rand(rng, Uniform(0, u.max_shift))
    theta = rand(rng, Uniform(0, π))
    phi = rand(rng, Uniform(0, 2π))
    return (
        radius * sin(theta) * cos(phi),
        radius * sin(theta) * sin(phi),
        radius * cos(theta)
    )
end
function Base.rand(rng::AbstractRNG, u::UniformBallShift{D, T}) where {D, T}
    error("Unsupported dimension $D for UniformBallShift")
end

##
abstract type ShiftMethod end
struct NoShift <: ShiftMethod end
marginal_shift(pp::PointSet, ::NoShift) = pp

## Toroidal shift
struct ToroidalShift{R <: Box, S <: Union{<:SpatialShift, <:NTuple}} <: ShiftMethod
    region::R
    shift::S
end
function ToroidalShift(box::Box)
    centered_box = inverse(Translate(to(centroid(box))...))(box)
    ToroidalShift(
        box,
        UniformShift(unitless_coords(centered_box.min), unitless_coords(centered_box.max))
    )
end
function Base.rand(rng::AbstractRNG, shift::ToroidalShift)
    ToroidalShift(shift.region, rand(rng, shift.shift))
end

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
function marginal_shift(pp::PointSet, shift_method::ToroidalShift)
    toroidal_shift(pp, shift_method.region, shift_method.shift)
end

## Standard shift
struct StandardShift{S <: Union{<:SpatialShift, <:NTuple}} <: ShiftMethod
    shift::S
end
"""
    StandardShift(shift)

A shift method that just applies a shift to all points. Does not make any assumptions about the region they are recorded on.
"""
StandardShift

Base.rand(rng::AbstractRNG, shift::StandardShift) = StandardShift(rand(rng, shift.shift))

function marginal_shift(pp::PointSet, shift_method::StandardShift)
    Translate(shift_method.shift...)(pp)
end

function marginal_shift(pp::PointPattern, shift_method)
    spatial_data(marginal_shift(observations(pp), shift_method), getregion(pp))
end

function shift_resample(
        data::MultipleSpatialDataTuple{P},
        statistic,
        shift_method::ShiftMethod,
        groups = 1:P;
        kwargs...
) where {P}
    shift_resample(
        Random.default_rng(),
        data,
        statistic,
        shift_method,
        groups;
        kwargs...
    )
end

function shift_resample(
        rng::AbstractRNG,
        data::MultipleSpatialDataTuple{P},
        statistic,
        shift_method::ShiftMethod,
        groups = 1:P;
        kwargs...
) where {P}
    @assert sort(reduce(vcat, groups))==1:P "groups of shifts should partition the space"
    group_shifts = Dict{eltype(groups), ShiftMethod}(
        group => rand(rng, shift_method) for group in groups
    )
    shifted_processes = spatial_data(
        ntuple(
            p -> observations(
                marginal_shift(data[p], group_shifts[findgroup(p, groups)])),
            Val{P}()),
        getregion(data))
    statistic(shifted_processes; kwargs...)
end

findgroup(p, groups) = groups[findfirst(g -> p ∈ g, groups)]
