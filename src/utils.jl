unitless_coords(p::Point) = getfield.(Meshes.CoordRefSystems._coords(coords(p)), :val)
unitless_spacing(grid::CartesianGrid) = getfield.(spacing(grid), :val)
unitless_origin(grid::CartesianGrid) = unitless_coords(grid.origin)
unitless_minimum(grid::CartesianGrid) = unitless_coords(minimum(grid))
unitless_maximum(grid::CartesianGrid) = unitless_coords(maximum(grid))
unitless_measure(region) = measure(region).val

"""
	points2coords(points::PointSet)

Converts PointSet to a tuple of view of reinterpreted arrays, whose `j'th entry is the `j`th coordinate of each of the points.
"""
function points2coords(points::PointSet)
    points2coords(parent(points), Val{typeof(to(points[1])[1].val)}(), Val{embeddim(points)}())
end
function points2coords(points, ::Val{T}, ::Val{D}) where {T, D}
    points_re = reinterpret(reshape, T, points)
	return ntuple(d -> view(points_re, d, :), Val{D}())
end
function points2coords(points, ::Val{T}, ::Val{1}) where {T}
    points_re = reinterpret(reshape, T, points)
	return (points_re,)
end

"""
	box2sides(box::Box)

Converts a box to a tuple of tuples, whose `j'th entry is the `j`th coordinate of the minimum and maximum of the box.
"""
function box2sides(box::Box)
	return ntuple(
		d -> (unitless_coords(minimum(box))[d], unitless_coords(maximum(box))[d]),
		embeddim(box),
	)
end

"""
	grid2side(g::CartesianGrid)

Converts a grid to a tuple of the location of the centers for each side.
"""
function grid2side(grid::CartesianGrid)
	return ntuple(
		d -> unitless_origin(grid)[d] .+ ((0.5:size(grid)[d]) .* unitless_spacing(grid)[d]),
		embeddim(grid),
	)
end

"""
	side2grid(sides)

inverse of grid2side
"""
function side2grid(sides::NTuple)
	start = first.(sides) .- step.(sides) ./ 2
	stop = last.(sides) .+ step.(sides) ./ 2
	n = length.(sides)
	CartesianGrid(start, stop, dims = n)
end

#### new utils

"""
	pad(x::AbstractArray{T,D}, n::NTuple{D,Int}) where {T,D}
	pad(x::AbstractArray{T,D}, n::Int) where {T,D}

Pad `x` with zeros to size `n`. 
Note that `size(x) ≤ n` must hold, and the new array is of size `n` not `size(x).+n`.
If an integer is provided for `n`, it is interpreted as `ntuple(d->n, Val{D}())`.
"""
pad(x::AbstractArray{T, D}, n::Int) where {T, D} = pad(x, ntuple(d -> n, Val{D}()))
function pad(x::AbstractArray{T, D}, n::NTuple{D, Int}) where {T, D}
	@assert all(size(x) .≤ n) "size(x) must be smaller than n"
	y = zeros(T, n)
	y[CartesianIndices(x)] .= x
	return y
end

"""
    centerpad(x, i)

Pads `x` with zeros on all sides by `i` elements.
For example, if `x = 1:2` and `i = 1`, then the output is `[0, 1, 2, 0]`.
"""
function centerpad(x, i)
    out = zeros(size(x) .+ 2 .* i)
    out[range.(1 .+ i, size(out).-i)...] .= x
    out
end

"""
	downsample(x::AbstractArray{T,D}, spacing) where {D,T}

Down sample an `AbstractArray`.
Specify `spacing` as an `Int`, for same spacing in all dims.
Specify `spacing` as an `NTuple{D,Int}` for different spacing.
Specify `spacing=nothing` to just return `x`.
"""
downsample(x::AbstractArray{T, D}, spacing::Int) where {D, T} =
	downsample(x, ntuple(d -> spacing, Val{D}()))
downsample(x::AbstractArray, ::Nothing) = x
function downsample(x::AbstractArray{T, D}, spacing::NTuple{D, Int}) where {D, T}
	@assert all(0 .< spacing .< size(x))
	ind = ntuple(d -> 1:spacing[d]:size(x, d), Val{D}())
	view(x, ind...)
end

"""
	upsample(x::AbstractArray, grid::CartesianGrid, freq_downsample)

Upsamples an array using linear interpolation, checking the grid size against x.
Note this is specific to downsample as used here, and shouldn't be used for general upsampling.
"""
upsample(x::AbstractArray, ::CartesianGrid, ::Nothing) = x
function upsample(x::AbstractArray, grid::CartesianGrid, freq_downsample)
	ndims(x) == embeddim(grid) ||
		throw(DimensionMismatch("x must have the same number of dimensions as the grid"))

	sides = grid2side(grid)
	@assert size(x) .* freq_downsample == length.(sides)
	x_sides = downsample.(sides, freq_downsample)
	x_intp = linear_interpolation(x_sides, x, extrapolation_bc = Interpolations.Line())
	return [x_intp(a...) for a in Iterators.product(sides...)]
end
