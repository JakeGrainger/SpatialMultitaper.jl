"""
    SpatialData{T,R}

A struct to hold spatial data and its associated region.

All methods in SpatialMultitaper.jl use `SpatialData` objects to represent spatial
datasets. Most will allow the user to pass other raw formats directly and then convert
them internally to `SpatialData` objects.

# Fields
- `data::T`: The spatial observations (PointSet, GeoTable, or collections thereof)
- `region::R`: The spatial region where data is observed

# Usage
You should not use this constructor directly.
You should specify data and region using the `spatial_data` function, which will
automatically unpack, validate, and mask the input data to the specified region.
"""
struct SpatialData{T, R <: Meshes.Geometry, N}
    data::T
    region::R
    names::N
end
Base.getindex(sd::SpatialData, idx) = SpatialData(sd.data[idx], sd.region, sd.names[idx])
Meshes.embeddim(sd::SpatialData) = embeddim(sd.region)
getregion(sd::SpatialData) = sd.region

"""
    spatial_data(data, region::R)

Creates a `SpatialData` object from raw spatial data and a specified region.
Automatically unpacks, validates, and masks the input data to the specified region.
"""
function spatial_data(data, region)
    formatted_data = unpack_observations(data)
    check_observations(formatted_data, region)
    masked_data = mask(formatted_data, region) # note this will copy
    names = defaultnames(masked_data)
    return SpatialData(masked_data, region, names)
end

ncol(sd::SpatialData{<:Tuple}) = length(sd.data)
ncol(sd::SpatialData{<:Vector}) = length(sd.data)
ncol(sd::SpatialData) = 1

"""
    observations(sd::SpatialData)

Return the spatial observations from a SpatialData object.
"""
observations(sd::SpatialData) = sd.data

Base.propertynames(sd::SpatialData) = sd.names

# type aliases to make it easy to refer to GeoTable types
const PointGeoTable = GeoTable{<:PointSet}
const GridGeoTable = GeoTable{<:CartesianGrid}

# The type aliases below represent every form of supported data input
const PointPattern = SpatialData{<:PointSet}
const MarkedPointPattern = SpatialData{<:PointGeoTable}
const GriddedData = SpatialData{<:GridGeoTable}
const MultipleSpatialDataTuple{P} = SpatialData{<:NTuple{P, Any}}
const MultipleSpatialDataVec = SpatialData{<:Vector}

# eventually will have proper names support, but for now downstream will always support names
defaultnames(data::Union{Tuple, Vector}) = 1:length(data)
defaultnames(data) = 1

# Helper function to ensure observations are in a supported format
abstract type SpatialDataUnpackMethod end
struct SingleUnpack <: SpatialDataUnpackMethod end
struct TupleUnpack <: SpatialDataUnpackMethod end
struct VectorUnpack <: SpatialDataUnpackMethod end

"""
    unpack_observations(data)

Unpacks input data into a supported format for spatial analysis.

# Behavior
- `PointSet` → returned as-is
- `GeoTable` with single column → returned as-is
- `GeoTable` with multiple columns → error for marked points, unpacked for grids
- `Tuple` of spatial data → recursively unpacked and flattened
- `Vector` of spatial data → recursively unpacked and concatenated

# Examples
```julia
# Single process
unpack_observations(pointset) # → pointset

# Multiple grid fields
unpack_observations(grid_geotable) # → (field1_geotable, field2_geotable, ...)

# Mixed collection
unpack_observations((pointset, grid)) # → (pointset, field1, field2, ...)
```

# Notes
Marked point processes with multiple marks are not supported - pass as separate
GeoTables instead.
"""
function unpack_observations(data::Union{PointSet, GeoTable})
    _unpack_observations(data, SingleUnpack())
end
function unpack_observations(data::NTuple{N, Any}) where {N}
    (_unpack_observations(data[1], TupleUnpack())...,
        unpack_observations(data[2:end])...)
end
unpack_observations(data::NTuple{1}) = _unpack_observations(data[1], TupleUnpack())
function unpack_observations(data::Vector)
    mapreduce(x -> _unpack_observations(x, VectorUnpack()), vcat, data)
end

_unpack_observations(data::PointSet, unpacktype) = _unpack_single(data, unpacktype)

function _unpack_observations(data::PointGeoTable, unpacktype)
    if ncol(data) - 1 > 1
        throw(ArgumentError(
            "Multi-mark point processes are not supported. Pass separate GeoTables " *
            "instead, though this is not recommended as processes should be jointly simple."
        ))
    end
    return _unpack_single(data, unpacktype)
end

function _unpack_observations(data::GridGeoTable, unpacktype)
    if ncol(data) - 1 > 1 # ncol includes the domain
        if unpacktype isa TupleUnpack()
            throw(ArgumentError(
                "Gridded datasets with >10 fields require Vector format. " *
                "Use Vector input or pass astuple=false."
            ))
        end
        if unpacktype isa SingleUnpack()
            return _unpack_multi(data, VectorUnpack())
        end
        return _unpack_multi(data, unpacktype)
    else
        return _unpack_single(data, unpacktype)
    end
end

_unpack_single(data, ::SingleUnpack) = data
_unpack_single(data, ::TupleUnpack) = (data,)
_unpack_single(data, ::VectorUnpack) = [data]
function _unpack_multi(data, ::TupleUnpack)
    tuple(single_geotable.(Ref(data), tuple(1:ncol(data)...))...)
end
_unpack_multi(data, ::VectorUnpack) = [single_geotable(data, i) for i in 1:ncol(data)]

"""
    single_geotable(data::GeoTable, idx::Integer)

Extract the `idx`-th column of a GeoTable as a new single-column GeoTable.
"""
function single_geotable(data::GeoTable, idx::Integer)
    name = propertynames(data)[idx]
    return georef(NamedTuple(name => getproperty(data, name),), domain(data))
end

"""
    check_observations(data, region)

Validate that spatial data is compatible with the specified region.

Ensures all data has the same parameter dimension as the region.
"""
check_observations(data::PointSet, region) = @argcheck embeddim(data) == embeddim(region)
function check_observations(data::GeoTable, region)
    @argcheck embeddim(domain(data)) == embeddim(region)
    @argcheck all(isfinite, values(data)[1])
end
function check_observations(data::Union{Tuple, Vector}, region)
    foreach(d -> check_observations(d, region), data)
end

"""
    mask(pp::PointSet, region)
    mask(geotable::GeoTable, region)

Masks a PointSet or GeoTable to only include points within the specified region.
If the data is already fully contained within the region, a copy is still made and returned.
If the data is a `GeoTable` on a grid, values outside the region are set to `NaN`.
"""
function mask(data::Union{Tuple, Vector}, region)
    map(d -> mask(d, region), data)
end

function mask(pp::PointSet, region)
    return PointSet([p for p in pp if p ∈ region])
end
function mask(geotable::PointGeoTable, region)
    pts = domain(geotable)
    marks = values(geotable)[1]
    name = propertynames(geotable)[1]

    newmarks = [m for (m, p) in zip(marks, pts) if p ∈ region]
    newpoints = mask(pts, region)
    return georef(NamedTuple((name => newmarks,)), newpoints)
end
function mask(geotable::GridGeoTable, region)
    data = values(geotable)
    grid = domain(geotable)

    newdata = deepcopy(data)
    for i in eachindex(newdata)
        for j in eachindex(newdata[i])
            if centroid(grid, j) ∉ region
                newdata[i][j] = NaN
            end
        end
    end
    return georef(newdata, grid)
end
