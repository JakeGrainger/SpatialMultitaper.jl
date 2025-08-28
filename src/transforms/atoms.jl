atom_estimate(data::PointSet, region) = length(data) / unitless_measure(region)
atom_estimate(data::GeoTable, region) = atom_estimate(domain(data), values(data)[1], region)
atom_estimate(domain::CartesianGrid, rf, region) = 0.0
atom_estimate(domain::PointSet, marks, region) =
    error("Atom for the marked case is not yet implemented")

atom_estimate(data::Union{Tuple,AbstractVector}, region) = atom_estimate.(data, Ref(region))
atom_estimate(data, region, ::UsualPartial) = atom_estimate(data, region)
atom_estimate(data, region, ::SplitPartial) = map(x -> 0.0, data) # split partial should have no atoms even on diagonals

# TODO: should link this to the existing mean methods