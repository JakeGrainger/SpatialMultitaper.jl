abstract type MeanEstimationMethod end

struct DefaultMean <: MeanEstimationMethod end

struct KnownMean{T} <: MeanEstimationMethod
    value::T
end

mean_estimate(data::PointSet, region, mean_method) = mean_estimate(
    data,
    [1.0 * (data[i] ∈ region) for i ∈ eachindex(data)],
    region,
    mean_method,
)
mean_estimate(data::GeoTable, region, mean_method) =
    mean_estimate(domain(data), values(data)[1], region, mean_method)

function mean_estimate(points::PointSet, marks, region, mean_method::DefaultMean)
    sum(marks[i] for i in eachindex(marks) if points[i] ∈ region) / unitless_measure(region)
end

function mean_estimate(grid::CartesianGrid, rf, region, mean_method::DefaultMean)
    mean(rf[i] for i in eachindex(rf) if centroid(grid, i) ∈ region)
end

function mean_estimate(domain, data, region, mean_method::KnownMean)
    mean_method.value
end
