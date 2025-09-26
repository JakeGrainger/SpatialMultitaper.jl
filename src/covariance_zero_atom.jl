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
covariance_zero_atom(data::NTuple{1}, region) = covariance_zero_atom(data[1], region)

function covariance_zero_atom(data::Tuple, region)
    diagm(SVector(covariance_zero_atom.(data, Ref(region))))
end

covariance_zero_atom(data::Vector, region) = diagm(covariance_zero_atom.(data, Ref(region)))

covariance_zero_atom(data::PointSet, region) = length(data) / unitless_measure(region)

function covariance_zero_atom(data::GeoTable, region)
    _covariance_zero_atom(domain(data), values(data)[1], region)
end

_covariance_zero_atom(::CartesianGrid, rf, region) = zero(eltype(rf))

_covariance_zero_atom(::PointSet, mark, region) = sum(abs2, mark) / unitless_measure(region)

# TODO: should link this to the existing mean methods
