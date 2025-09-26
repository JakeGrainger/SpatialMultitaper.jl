abstract type MeanEstimationMethod end

struct DefaultMean <: MeanEstimationMethod end
Base.getindex(m::DefaultMean, i) = DefaultMean()
Base.lastindex(m::DefaultMean) = 1

struct KnownMean{T} <: MeanEstimationMethod
    value::T
end
Base.getindex(m::KnownMean, i) = KnownMean(m.value[i])
Base.lastindex(m::KnownMean) = lastindex(m.value)

checkmeanmethod(data, ::DefaultMean) = nothing
function checkmeanmethod(data, m::KnownMean)
    if length(data) != length(m.value)
        throw(ArgumentError("data and mean_method must have the same length, got $(length(data)) and $(length(m.value))."))
    end
end

function mean_estimate(data::NTuple{D}, region, mean_method) where {D}
    SVector(ntuple(j -> _mean_estimate(data[j], region, mean_method[j]), Val{D}()))
end

function mean_estimate(data::AbstractVector, region, mean_method)
    [mean_estimate(data[i], region, mean_method[i]) for i in eachindex(data)]
end

function mean_estimate(data::NTuple{1}, region, mean_method)
    _mean_estimate(data[1], region, mean_method[1])
end

function mean_estimate(data::Union{PointSet, GeoTable}, region, mean_method)
    _mean_estimate(data, region, mean_method)
end

##
function _mean_estimate(data, region, mean_method::KnownMean)
    mean_method.value
end

function _mean_estimate(data::PointSet, region, ::DefaultMean)
    sum(x -> x ∈ region, data) / unitless_measure(region)
end
function _mean_estimate(data::GeoTable, region, ::DefaultMean)
    _mean_estimate(domain(data), values(data)[1], region, DefaultMean())
end

function _mean_estimate(points::PointSet, marks::AbstractVector, region, ::DefaultMean)
    sum(marks[i] for i in eachindex(marks) if points[i] ∈ region) / unitless_measure(region)
end

function _mean_estimate(grid::CartesianGrid, rf, region, ::DefaultMean)
    mean(rf[i] for i in eachindex(rf) if centroid(grid, i) ∈ region)
end
