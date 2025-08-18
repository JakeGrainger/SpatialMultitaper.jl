atom_estimate(data::PointSet, region) = length(data) / unitless_measure(region)
atom_estimate(data::GeoTable, region) = atom_estimate(domain(data), values(data)[1], region)
atom_estimate(domain::CartesianGrid, rf, region) = 0.0
atom_estimate(domain::PointSet, marks, region) =
    error("Atom for the marked case is not yet implemented")

atom_estimate(data::Union{Tuple,AbstractVector}, region) = atom_estimate.(data, Ref(region))
# TODO: should link this to the existing mean methods