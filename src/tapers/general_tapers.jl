"""
	optimaltapers(region::Meshes.GeometryOrDomain, grid::CartesianGrid; wavenumber_region::Meshes.GeometryOrDomain, ntapers::Int, wavenumber_res, wavenumber_downsample=nothing, real_tapers=true, tol=0.0)

Function to compute the optimal tapers for a given `region`, `grid` and region in wavenumber.

# Arguments:
- `region`: The spatial observation region. Should be of type `Geometry`.
- `grid`: A `CartesianGrid` containing the region on which we wish to have a taper.
- `wavenumber_region`: The region in wavenumber we want to concentrate the taper on. Typically a `Ball` centered at zero.
- `ntapers`: The number of desired tapers.
- `wavenumber_res`: The oversampling to be used in wavenumber.
- `real_tapers`: Optional argument to decided if real or complex tapers should be provided. Default is `true` which provides real tapers.
- `tol`: Optional argument passed to the `eigs` function of `Arpack`. You likely need to play with this.
- `check_grid`: Optional argument to decide if the grid should be checked for compatibility with the region. Default is `true`. Sometimes due to float error this can fail, so you may want to set this to `false` if you are sure the grid and region are compatible.

# Note about errors:
If you have an error, it is likely a convergence problem. Try setting a larger value for `tol`.
"""
function optimaltapers(
        region::Meshes.Geometry,
        grid::CartesianGrid;
        wavenumber_region::Meshes.Geometry,
        ntapers::Int,
        wavenumber_res,
        wavenumber_downsample = nothing,
        real_tapers = true,
        tol = 0.0,
        check_grid = true
)
    @assert embeddim(grid)==embeddim(region) "The region and grid must have the same number of dimensions."
    if check_grid
        @assert boundingbox(region)⊆boundingbox(grid) "The region of interest is not covered by the grid, please provide a bigger grid or smaller region. Note this could be due to floating point error."
    end
    @assert embeddim(wavenumber_region)==embeddim(region) "The region and wavenumber_region must have the same number of dimensions."
    wavenumber_res = checkwavenumberres(grid, wavenumber_res)

    R = padto(
        downsample(pixelate_region(grid, region), wavenumber_downsample), wavenumber_res)
    K = pixelate_region(
        fftfreq.(wavenumber_res, 1 ./ downsample_spacing(grid, wavenumber_downsample)),
        wavenumber_region
    )
    h_oversize, λ = compute_eigenfunctions(
        R, K, ntapers; real_tapers = real_tapers, tol = tol)
    h = [prod(sqrt, unitless_spacing(grid)) .*
         reprocess(h_oversize[i], grid, wavenumber_downsample)
         for i in eachindex(h_oversize)]
    return h, λ
end

function checkwavenumberres(grid::CartesianGrid, wavenumber_res)
    if !(wavenumber_res isa Int)
        @assert length(wavenumber_res)==embeddim(grid) "wavenumber_res must be an integer or a tuple of integers with the same length as the dimension of grid."
    end
    if all(size(grid) .≤ wavenumber_res)
        return wavenumber_res
    else
        @warn "wavenumber_res is smaller than grid size, will set to the grid size but you may wish to increase this."
        return max.(size(grid), wavenumber_res)
    end
end

"""
	reprocess(h_large, grid, wavenumber_downsample)

Undoes the downsampling and padding.
"""
function reprocess(h_large::AbstractArray{T, D}, grid, wavenumber_downsample) where {D, T}
    return newweight(h_large, wavenumber_downsample) .* upsample(
        h_large[ntuple(
            d -> 1:downsample_size(grid, wavenumber_downsample)[d], Val{D}())...],
        grid,
        wavenumber_downsample
    )
end

"""
	newweight(h, wavenumber_downsampling)

Necessary because eigs gives normalised vector under the downsampling.
"""
newweight(::AbstractArray, ::Nothing) = 1.0
function newweight(::AbstractArray{T, D}, wavenumber_downsample::Real) where {T, D}
    sqrt(1 / wavenumber_downsample)^D
end
function newweight(
        ::AbstractArray{T, D}, wavenumber_downsample::NTuple{D, Real}) where {T, D}
    prod(sqrt, inv.(wavenumber_downsample))
end

downsample_size(grid::CartesianGrid, ::Nothing) = size(grid)
function downsample_size(grid::CartesianGrid, wavenumber_downsample)
    size(grid) .÷ wavenumber_downsample
end

downsample_spacing(grid::CartesianGrid, ::Nothing) = unitless_spacing(grid)
function downsample_spacing(grid::CartesianGrid, wavenumber_downsample)
    unitless_spacing(grid) .* wavenumber_downsample
end

function pixelate_region(grid::CartesianGrid, shp::Meshes.Geometry)
    Δ = spacing(grid)
    b = Box(Point(.-Δ), Point(Δ))
    R = [Translate(unitless_coords(centroid(grid, i)))(b) ⊆ shp for i in eachindex(grid)]
    return reshape(R, size(grid))
end

# when doing this in wavenumber
function pixelate_region(g::NTuple{D, Frequencies}, shp::Meshes.GeometryOrDomain) where {D}
    R = [Point(i) ∈ shp for i in Iterators.product(g...)]
    return R
end
