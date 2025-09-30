"""
    unitless_coords(p::Point)

Return the coordinates of a point without units.
"""
unitless_coords(p::Point) = Meshes.ustrip.(Meshes.CoordRefSystems._coords(coords(p)))

"""
    unitless_spacing(grid::CartesianGrid)

Return the spacing of a CartesianGrid without units.
"""
unitless_spacing(grid::CartesianGrid) = Meshes.ustrip.(spacing(grid))

"""
    unitless_origin(grid::CartesianGrid)

Return the origin of a CartesianGrid without units.
"""
unitless_origin(grid::CartesianGrid) = unitless_coords(grid.origin)

"""
    unitless_minimum(grid::CartesianGrid)

Return the minimum coordinates of a CartesianGrid without units.
"""
unitless_minimum(grid::CartesianGrid) = unitless_coords(minimum(grid))

"""
    unitless_maximum(grid::CartesianGrid)

Return the maximum coordinates of a CartesianGrid without units.
"""
unitless_maximum(grid::CartesianGrid) = unitless_coords(maximum(grid))

"""
    unitless_measure(region)

Return the measure of a `region` without units.
"""
unitless_measure(region) = Meshes.ustrip(measure(region))

"""
    point2unitlesstype(p::Point)

Return the type of the unitless coordinates of a Point.
"""
point2unitlesstype(p::Point) = typeof(to(p)[1].val)

"""
	points2coords(points::PointSet)

Return a `Tuple`, with `j'th entry a vector of the `j`th entry of each `point` in `points`.
"""
function points2coords(points::PointSet)
    _points2coords(
        parent(points),
        Val{point2unitlesstype(points[1])}(),
        Val{embeddim(points)}()
    )
end
function _points2coords(points, ::Val{T}, ::Val{D}) where {T, D}
    points_re = reinterpret(reshape, T, points)
    return ntuple(d -> view(points_re, d, :), Val{D}())
end
function _points2coords(points, ::Val{T}, ::Val{1}) where {T}
    points_re = reinterpret(reshape, T, points)
    return (points_re,)
end

"""
	box2sides(box::Box)

Converts a box to a `Tuple` of `Tuple`s, whose `j'th entry is min and max of the `j`th side.
"""
function box2sides(box::Box)
    return ntuple(
        d -> (unitless_coords(minimum(box))[d], unitless_coords(maximum(box))[d]),
        embeddim(box)
    )
end

"""
	grid2sides(g::CartesianGrid)

Converts a grid to a tuple of the location of the centers for each side.
"""
function grid2sides(grid::CartesianGrid)
    vertex_sides = Meshes.xyz(grid)
    return map(_grid2side, vertex_sides)
end
function _grid2side(side_with_units)
    side = range(Meshes.ustrip(first(side_with_units)),
        Meshes.ustrip(last(side_with_units)), step = Meshes.ustrip(step(side_with_units)))

    @assert length(side) == length(side_with_units)
    (side .+ step(side) / 2)[1:(end - 1)]
end

# function grid2sides(grid::CartesianGrid, ::Val{T}, ::Val{D}) where {T, D}
#     return ntuple(
#         d -> unitless_origin(grid)[d] .+
#              ((convert(T, 0.5):size(grid)[d]) .* unitless_spacing(grid)[d]),
#         Val{D}()
#     )
# end
# """
# 	sides2grid(sides)

# inverse of grid2sides
# """
# function sides2grid(sides::NTuple)
#     start = first.(sides) .- step.(sides) ./ 2
#     stop = last.(sides) .+ step.(sides) ./ 2
#     n = length.(sides)
#     CartesianGrid(start, stop, dims = n)
# end

@deprecate sides2grid(sides::NTuple) nothing

"""
	padto(x::AbstractArray{T,D}, n::NTuple{D,Int}) where {T,D}
	padto(x::AbstractArray{T,D}, n::Int) where {T,D}

Pad `x` with zeros to size `n`.
Important: the resultant array is size `n`, not `size(x) .+ n`.
Note that `size(x) ≤ n` must hold, and the new array is of size `n` not `size(x).+n`.
If an integer is provided for `n`, it is interpreted as the same size in all dimensions.
"""
padto(x::AbstractArray{T, D}, n::Int) where {T, D} = padto(x, ntuple(Returns(n), Val{D}()))
function padto(x::AbstractArray{T, D}, n::NTuple{D, Int}) where {T, D}
    if !all(size(x) .≤ n)
        err = ArgumentError("""size(x) must be smaller than n in each dimension,
                                but size(x)=$(size(x)) and n=$(n)"""
        )
        throw(err)
    end
    y = zeros(T, n)
    y[CartesianIndices(x)] .= x
    return y
end

"""
    centerpad(x::AbstractArray, i)

Pads `x` with zeros on all sides by `i` elements.
For example, if `x = 1:2` and `i = 1`, then the output is `[0, 1, 2, 0]`.
"""
function centerpad(x::AbstractArray, i)
    if !all(i .>= 0)
        throw(ArgumentError("i must be non-negative, but got i=$(i)"))
    end
    if !all(isa.(i, Int))
        throw(ArgumentError("i must be an integer or tuple of integers, but got i=$(i)"))
    end
    if length(i) > 1 && length(i) != ndims(x)
        err = ArgumentError(
            "i must be a single integer or have length ndims(x)=$(ndims(x)), but got i=$(i)"
        )
        throw(err)
    end

    out = zeros(eltype(x), size(x) .+ 2 .* i)
    out[range.(1 .+ i, size(out) .- i)...] .= x
    out
end

"""
	downsample(x::AbstractArray{T,D}, spacing) where {D,T}

Down sample an `AbstractArray`.
Specify `spacing` as an `Int`, for same spacing in all dims.
Specify `spacing` as an `NTuple{D,Int}` for different spacing.
Specify `spacing=nothing` to just return `x`.
"""
function downsample(x::AbstractArray{T, D}, spacing::Int) where {D, T}
    downsample(x, ntuple(Returns(spacing), Val{D}()))
end
downsample(x::AbstractArray, ::Nothing) = x
function downsample(x::AbstractArray{T, D}, spacing::NTuple{D, Int}) where {D, T}
    if !all(0 .< spacing .< size(x))
        err = ArgumentError(
            "spacing must be positive and less than size(x)=$(size(x)), but got spacing=$(spacing)"
        )
        throw(err)
    end
    ind = ntuple(d -> firstindex(x, d):spacing[d]:lastindex(x, d), Val{D}())
    view(x, ind...)
end

"""
	upsample(x::AbstractArray, grid::CartesianGrid, wavenumber_downsample)

Upsamples an array using linear interpolation, checking the grid size against x.
Note this is specific to downsample as used here, and shouldn't be used for general upsampling.
"""
upsample(x::AbstractArray, ::CartesianGrid, ::Nothing) = x
function upsample(x::AbstractArray, grid::CartesianGrid, wavenumber_downsample)
    if ndims(x) != embeddim(grid)
        throw(DimensionMismatch("x must have the same number of dimensions as the grid"))
    end

    sides = grid2sides(grid)
    if size(x) .* wavenumber_downsample != length.(sides)
        err = DimensionMismatch(
            """we need size(x) .* wavenumber_downsample == length.(sides), but size(x)=$(size(x)),
            wavenumber_downsample=$(wavenumber_downsample), and length.(sides)=$(length.(sides))"""
        )
        throw(err)
    end
    x_sides = downsample.(sides, wavenumber_downsample)
    x_intp = linear_interpolation(x_sides, x, extrapolation_bc = Interpolations.Line())
    return [x_intp(a...) for a in Iterators.product(sides...)]
end
