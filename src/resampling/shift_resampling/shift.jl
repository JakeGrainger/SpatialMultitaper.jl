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
ToroidalShift(data::SpatialData) = ToroidalShift(getregion(data))
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

findgroup(p, groups) = groups[findfirst(g -> p ∈ g, groups)]

function apply_shifts(data::MultipleSpatialDataTuple{P}, shifts, groups) where {P}
    return spatial_data(
        ntuple(p -> observations(marginal_shift(data[p], shifts[findgroup(p, groups)])),
            Val{P}()),
        getregion(data))
end

function apply_shifts(data::MultipleSpatialDataVec, shifts, groups)
    return spatial_data(
        [observations(marginal_shift(data[p], shifts[findgroup(p, groups)]))
         for p in 1:ncols(data)],
        getregion(data))
end

function apply_shifts(data::SingleProcessData, shifts, groups)
    return spatial_data(
        marginal_shift(observations(data), shifts[groups[1]]), getregion(data))
end

##
struct ShiftResampler{T, S, R, G, M, K}
    data::R
    shift_method::S
    groups::G
    mem::M
    kwargs::K
    function ShiftResampler{T}(
            data::R, shift_method::S, groups::G, mem::M, kwargs::K) where {T, S, R, G, M, K}
        @argcheck sort(reduce(vcat, groups)) == 1:ncol(data)
        new{T, S, R, G, M, K}(data, shift_method, groups, mem, kwargs)
    end
end

function ShiftResampler(statistic, data::SpatialData,
        shift_method::ShiftMethod = default_shift_method(data),
        groups = 1:ncol(data); kwargs...)
    T = functional_statistic_type(statistic, data)
    return ShiftResampler(T, data, shift_method, groups; kwargs...)
end

function ShiftResampler(::Type{T}, data::SpatialData,
        shift_method::ShiftMethod = default_shift_method(data),
        groups = 1:ncol(data); kwargs...) where {T}
    resolved_kwargs = resolve_parameters(T, data; kwargs...)
    mem = preallocate_memory(T, data; resolved_kwargs...)
    ShiftResampler{T}(data, shift_method, groups, mem, resolved_kwargs)
end

default_shift_method(data::SpatialData) = default_shift_method(getregion(data))
default_shift_method(region::Box) = ToroidalShift(region)

function shift_resample!(rng::AbstractRNG, resampler::ShiftResampler{T}) where {T}
    groups = resampler.groups
    group_shifts = Dict(group => rand(rng, resampler.shift_method) for group in groups)

    shifted_data = apply_shifts(resampler.data, group_shifts, groups)
    result = estimate_function!(T, resampler.mem, shifted_data; resampler.kwargs...)
    return deepcopy(result) # ensure no references to mem
end

function shift_resample!(resampler::ShiftResampler)
    return shift_resample!(Random.default_rng(), resampler)
end
