"""
    covariance_zero_atom(data, region)

Estimate the zero atom of the reduced covariance measure for spatial data.

For point data, computes the density as the number of points divided by the region's
measure. For marked point data, computes the sum of squared mark values divided by
the region's measure. For gridded data, returns zero.

Returns a scalar for single datasets or a diagonal matrix for multiple datasets.

# Arguments
- `data`: Spatial data (PointSet, GeoTable, CartesianGrid, or collections thereof)
- `region`: Spatial region for which to compute the estimate

# Returns
- Scalar estimate for single datasets
- Diagonal matrix for multiple datasets (Tuple or Vector of datasets)
"""
function covariance_zero_atom(data::MultipleSpatialDataTuple{P}) where {P}
    diagm(SVector(ntuple(i -> covariance_zero_atom(data[i]), Val{P}())))
end

function covariance_zero_atom(data::MultipleSpatialDataVec)
    diagm([covariance_zero_atom(data[i]) for i in 1:ncol(data)])
end

function covariance_zero_atom(data::PointPattern)
    length(observations(data)) / unitless_measure(getregion(data))
end

function covariance_zero_atom(data::GriddedData)
    rf = values(observations(data))[1]
    zero(eltype(rf))
end

function covariance_zero_atom(data::MarkedPointPattern)
    mark = values(observations(data))[1]
    sum(abs2, mark) / unitless_measure(getregion(data))
end

# TODO: should link this to the existing mean methods
