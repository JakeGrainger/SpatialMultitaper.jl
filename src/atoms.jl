atom_estimate(data::PointSet, region) = length(data) / unitless_measure(region)
atom_estimate(data::GeoTable, region) = atom_estimate(domain(data), values(data)[1], region)
atom_estimate(domain::CartesianGrid, rf, region) = zero(eltype(rf))
function atom_estimate(domain::PointSet, marks::AbstractVector, region)
    sum(abs2, marks) / unitless_measure(region)
end

atom_estimate(data::NTuple{1}, region) = atom_estimate(data[1], region)
atom_estimate(data::Tuple, region) = diagm(SVector(atom_estimate.(data, Ref(region))))
atom_estimate(data::Vector, region) = diagm(atom_estimate.(data, Ref(region)))
atom_estimate(data, region) = atom_estimate(data, region)

# TODO: should link this to the existing mean methods
