module TestData

using SpatialMultitaper, GeoStatsProcesses, StableRNGs

import Random: AbstractRNG

# Helper functions for arbitrary dimensions
_default_region_min(ndims) = ntuple(_ -> -1.0, ndims)
function _default_region_max(ndims)
    ndims == 1 ? (4.0,) : ndims == 2 ? (3.0, 5.0) : ntuple(i -> i + 2.0, ndims)
end
function _default_grid_dims(ndims)
    ndims == 1 ? (25,) : ndims == 2 ? (21, 14) : ntuple(i -> 20 - i, ndims)
end

function _make_region(region_min, region_max)
    return Box(Point(region_min), Point(region_max))
end

function _make_grid(rng, region, grid_dims)
    total_points = prod(grid_dims)
    return georef((rf = rand(rng, total_points),),
        CartesianGrid(minimum(region), maximum(region), dims = grid_dims))
end

"""
    make_points_example(rng::AbstractRNG; n_processes = 2, as_vector = false,
        dim = 2, region_min = _default_region_min(dim),
        region_max = _default_region_max(dim), point_number = 100)
                       region_min = _default_region_min(2),
                       region_max = _default_region_max(2),
                       point_number = 100)

Generate test data consisting of point processes (unmarked point patterns).

Creates `n_processes` independent Poisson point processes within a rectangular region,
useful for testing spatial multitaper methods on point pattern data in arbitrary dimensions.

# Arguments
- `rng::AbstractRNG`: Random number generator for reproducible results
- `n_processes::Int = 2`: Number of independent point processes to generate
- `as_vector::Bool = false`: Return format control
  - `false`: Returns `((pattern...,), region)` tuple format
  - `true`: Returns `[pattern...]` vector format
- `region_min`: Minimum coordinates as tuple (default: `(-1, -1)` for 2D)
- `region_max`: Maximum coordinates as tuple (default: `(3, 5)` for 2D)
- `point_number::Real = 100`: Average number of points total

# Returns
- If `as_vector = false`: `((patterns...,), region)` where patterns is a tuple of `PointSet`s
- If `as_vector = true`: `[patterns...]` vector of `PointSet`s
- `region`: `Box` defining the spatial domain (only returned when `as_vector = false`)

# Examples
```julia
using StableRNGs
rng = StableRNG(123)

# Default 2D case
patterns, region = make_points_example(rng, n_processes=3)

# 1D case
patterns = make_points_example(rng, region_min=(-2.0,), region_max=(4.0,), as_vector=true)

# 3D case with custom density
patterns = make_points_example(rng, region_min=(-1,-1,-1), region_max=(2,2,2),
                              point_number=50, as_vector=true)
```

# Notes
- Supports arbitrary spatial dimensions (1D, 2D, 3D, etc.)
- Point density is per unit volume (length/area/volume/hypervolume)
- All processes are independent realizations of the same Poisson process
- Default region: `Box(Point(-1, -1), Point(3, 5))` for backward compatibility
"""
function make_points_example(rng::AbstractRNG; n_processes = 2, return_type = :tuple,
        dim = 2, region_min = _default_region_min(dim),
        region_max = _default_region_max(dim), point_number = 100)
    region = _make_region(region_min, region_max)
    pp_model = PoissonProcess(point_number / measure(region).val)
    pattern = rand(rng, pp_model, region, n_processes)

    if return_type == :vector
        return spatial_data(pattern, region)
    elseif return_type == :tuple
        return spatial_data((pattern...,), region)
    elseif return_type == :single
        @assert n_processes==1 "Single return type requires n_processes=1"
        return spatial_data(pattern[1], region)
    elseif return_type == :raw
        return pattern, region
    else
        throw(ArgumentError("Invalid return_type: $return_type"))
    end
end

"""
    make_grids_example(rng::AbstractRNG; n_processes = 2, return_type = :tuple,
                      dim = 2, region_min = _default_region_min(dim),
                      region_max = _default_region_max(dim),
                      grid_dims = _default_grid_dims(dim))

Generate test data consisting of gridded spatial data (regular lattice data).

Creates `n_processes` independent random fields on regular Cartesian grids within a
rectangular region, useful for testing spatial multitaper methods on gridded data
in arbitrary dimensions.

# Arguments
- `rng::AbstractRNG`: Random number generator for reproducible results
- `n_processes::Int = 2`: Number of independent gridded processes to generate
- `return_type::Symbol = :tuple`: Return format control
  - `:tuple`: Returns `SpatialData` with tuple format
  - `:vector`: Returns `SpatialData` with vector format
  - `:single`: Returns single `SpatialData` (requires `n_processes=1`)
  - `:raw`: Returns `(griddata, region)` without `SpatialData` wrapper
- `dim::Int = 2`: Spatial dimension (1D, 2D, 3D, etc.)
- `region_min`: Minimum coordinates as tuple (default based on `dim`)
- `region_max`: Maximum coordinates as tuple (default based on `dim`)
- `grid_dims`: Grid dimensions as tuple (default based on `dim`)

# Returns
Depends on `return_type` - see argument description above.

# Examples
```julia
using StableRNGs
rng = StableRNG(123)

# Default 2D case
grids = make_grids_example(rng, n_processes=2)

# 1D case with custom grid
grids = make_grids_example(rng, dim=1, grid_dims=(100,), return_type=:vector)

# 3D case
grids = make_grids_example(rng, dim=3, grid_dims=(20,20,20), return_type=:single, n_processes=1)
```

# Notes
- Supports arbitrary spatial dimensions (1D, 2D, 3D, etc.)
- Each grid contains independent random values drawn from uniform distribution
- Grid data is georeferenced using Meshes.jl/GeoStats.jl framework
"""
function make_grids_example(rng::AbstractRNG; n_processes = 2, return_type = :tuple,
        dim = 2, region_min = _default_region_min(dim),
        region_max = _default_region_max(dim),
        grid_dims = _default_grid_dims(dim))
    region = _make_region(region_min, region_max)
    griddata = [_make_grid(rng, region, grid_dims) for _ in 1:n_processes]

    if return_type == :vector
        return spatial_data(griddata, region)
    elseif return_type == :tuple
        return spatial_data((griddata...,), region)
    elseif return_type == :single
        @assert n_processes==1 "Single return type requires n_processes=1"
        return spatial_data(griddata[1], region)
    elseif return_type == :raw
        return griddata, region
    else
        throw(ArgumentError("Invalid return_type: $return_type"))
    end
end

"""
    make_marked_example(rng::AbstractRNG; n_processes = 2, return_type = :tuple,
                       dim = 2, region_min = _default_region_min(dim),
                       region_max = _default_region_max(dim),
                       point_number = 100)

Generate test data consisting of marked point processes (point patterns with associated data).

Creates `n_processes` independent marked point processes by first generating point patterns
and then attaching random marks (scalar values) to each point. Useful for testing spatial
multitaper methods on marked point pattern data in arbitrary dimensions.

# Arguments
- `rng::AbstractRNG`: Random number generator for reproducible results
- `n_processes::Int = 2`: Number of independent marked point processes to generate
- `return_type::Symbol = :tuple`: Return format control (see `make_points_example`)
- `dim::Int = 2`: Spatial dimension (1D, 2D, 3D, etc.)
- `region_min`: Minimum coordinates as tuple (default based on `dim`)
- `region_max`: Maximum coordinates as tuple (default based on `dim`)
- `point_number::Real = 100`: Average number of points total

# Returns
Depends on `return_type` - marked point processes with associated random mark values.

# Examples
```julia
using StableRNGs
rng = StableRNG(123)

# Default 2D case
marked_patterns = make_marked_example(rng, n_processes=3, return_type=:vector)

# 1D case
marked_1d = make_marked_example(rng, dim=1, return_type=:single, n_processes=1)
```

# Notes
- Built on top of `make_points_example` - supports arbitrary dimensions
- Marks are independent random values from uniform distribution [0,1]
- Each point gets exactly one scalar mark value
"""
function make_marked_example(rng::AbstractRNG; n_processes = 2, return_type = :tuple,
        dim = 2, region_min = _default_region_min(dim),
        region_max = _default_region_max(dim),
        point_number = 100)
    points, region = make_points_example(
        rng, n_processes = n_processes, return_type = :raw,
        dim = dim, region_min = region_min, region_max = region_max,
        point_number = point_number)
    marked_points = map(x -> georef((mark = rand(rng, length(x)),), x), points)

    if return_type == :vector
        return spatial_data(marked_points, region)
    elseif return_type == :tuple
        return spatial_data((marked_points...,), region)
    elseif return_type == :single
        @assert n_processes==1 "Single return type requires n_processes=1"
        return spatial_data(marked_points[1], region)
    elseif return_type == :raw
        return marked_points, region
    else
        throw(ArgumentError("Invalid return_type: $return_type"))
    end
end

"""
    make_mixed_example(rng::AbstractRNG; n_processes = (1, 1, 1), return_type = :tuple,
                      dim = 2, region_min = _default_region_min(dim),
                      region_max = _default_region_max(dim),
                      grid_dims = _default_grid_dims(dim),
                      point_number = 100)

Generate test data with mixed process types (points, grids, and marked points combined).

Creates a heterogeneous collection of spatial processes by combining point processes,
gridded data, and marked point processes. This is useful for testing multitaper methods
on mixed data types that might occur in real applications, in arbitrary dimensions.

# Arguments
- `rng::AbstractRNG`: Random number generator for reproducible results
- `n_processes::Tuple{Int,Int,Int} = (1, 1, 1)`: Number of each process type as `(points, grids, marks)`
- `return_type::Symbol = :tuple`: Return format control
- `dim::Int = 2`: Spatial dimension (1D, 2D, 3D, etc.)
- `region_min`: Minimum coordinates as tuple (default based on `dim`)
- `region_max`: Maximum coordinates as tuple (default based on `dim`)
- `grid_dims`: Grid dimensions as tuple (default based on `dim`)
- `point_number::Real = 100`: Average number of points total

# Returns
Depends on `return_type` - mixed collection of all process types sharing the same region.

# Examples
```julia
using StableRNGs
rng = StableRNG(123)

# Default 2D case
mixed_data = make_mixed_example(rng, n_processes=(2, 1, 3))

# 1D case with custom parameters
mixed_1d = make_mixed_example(rng, dim=1, n_processes=(2, 2, 2),
                             grid_dims=(100,), return_type=:vector)
```

# Notes
- Combines outputs from other make_*_example functions
- Supports arbitrary spatial dimensions
- Order in output: unmarked points, then grids, then marked points
"""
function make_mixed_example(
        rng::AbstractRNG; n_processes = (1, 1, 1), return_type = :tuple,
        dim = 2, region_min = _default_region_min(dim),
        region_max = _default_region_max(dim),
        grid_dims = _default_grid_dims(dim),
        point_number = 100)
    points, region = make_points_example(
        rng, n_processes = n_processes[1], return_type = :raw,
        dim = dim, region_min = region_min, region_max = region_max,
        point_number = point_number)
    grids, region_2 = make_grids_example(
        rng, n_processes = n_processes[2], return_type = :raw,
        dim = dim, region_min = region_min, region_max = region_max,
        grid_dims = grid_dims)
    marks, region_3 = make_marked_example(
        rng, n_processes = n_processes[3], return_type = :raw,
        dim = dim, region_min = region_min, region_max = region_max,
        point_number = point_number)
    @assert region==region_2==region_3 "Regions must match"

    all_processes = vcat(points, grids, marks)

    if return_type == :vector
        return spatial_data(all_processes, region)
    elseif return_type == :tuple
        return spatial_data((all_processes...,), region)
    elseif return_type == :single
        @assert length(all_processes)==1 "Single return type requires total n_processes=1"
        return spatial_data(all_processes[1], region)
    elseif return_type == :raw
        return all_processes, region
    else
        throw(ArgumentError("Invalid return_type: $return_type"))
    end
end

"""
    basic_test_grid(T, min, max, dims)

Create a test `CartesianGrid` with specified type `T`, min, max, and dims.

# Arguments
- `T`: The numeric type for coordinates
- `min`: Minimum coordinates as tuple
- `max`: Maximum coordinates as tuple
- `dims`: Grid dimensions as tuple

# Returns
A `CartesianGrid` with the specified parameters
"""
function basic_test_grid(T, min, max, dims)
    return CartesianGrid(convert.(T, min), convert.(T, max), dims = dims)
end

export make_points_example, make_grids_example, make_marked_example, make_mixed_example,
       basic_test_grid

end # module TestData
