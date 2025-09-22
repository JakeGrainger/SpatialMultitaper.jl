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

function mean_estimate(data::PointSet, region, mean_method)
    mean_estimate(
        data,
        [1.0 * (data[i] ∈ region) for i in eachindex(data)],
        region,
        mean_method
    )
end
function mean_estimate(data::GeoTable, region, mean_method)
    mean_estimate(domain(data), values(data)[1], region, mean_method)
end

function mean_estimate(points::PointSet, marks::AbstractVector, region, ::DefaultMean)
    sum(marks[i] for i in eachindex(marks) if points[i] ∈ region) / unitless_measure(region)
end

function mean_estimate(grid::CartesianGrid, rf, region, ::DefaultMean)
    mean(rf[i] for i in eachindex(rf) if centroid(grid, i) ∈ region)
end

function mean_estimate(data::AbstractVector, region, ::DefaultMean)
    mean_estimate.(data, Ref(region), Ref(DefaultMean()))
end

function mean_estimate(data::Tuple, region, ::DefaultMean)
    SVector(mean_estimate.(data, Ref(region), Ref(DefaultMean())))
end

function mean_estimate(domain, data, region, mean_method::KnownMean)
    mean_method.value
end
