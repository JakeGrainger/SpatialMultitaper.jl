atom_estimate(data::PointSet, region) = length(data) / unitless_measure(region)
atom_estimate(data::GeoTable, region) = atom_estimate(domain(data), values(data)[1], region)
atom_estimate(domain::CartesianGrid, rf, region) = zero(eltype(rf))
function atom_estimate(domain::PointSet, marks::AbstractVector, region)
    @warn("atom_estimate for the marked case is not yet implemented")
    return zero(eltype(marks))
end

atom_estimate(data::NTuple{1}, region) = atom_estimate(data[1], region)
atom_estimate(data::Tuple, region) = diagm(SVector(atom_estimate.(data, Ref(region))))
atom_estimate(data::Vector, region) = diagm(atom_estimate.(data, Ref(region)))
atom_estimate(data, region) = atom_estimate(data, region)

# TODO: should link this to the existing mean methods
