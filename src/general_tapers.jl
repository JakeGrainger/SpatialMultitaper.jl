"""
	optimaltapers(region::Meshes.GeometryOrDomain, grid::CartesianGrid; freq_region::Meshes.GeometryOrDomain, ntapers::Int, freq_res, freq_downsample=nothing, real_tapers=true, tol=0.0)

Function to compute the optimal tapers for a given `region`, `grid` and region in frequency.

# Arguments:
- `region`: The spatial observation region. Should be of type `Geometry` or `Domain`.
- `grid`: A `CartesianGrid` containing the region on which we wish to have a taper.
- `freq_region`: The region in frequency we want to concentrate the taper on. Typically a `Ball` centered at zero.
- `ntapers`: The number of desired tapers.
- `freq_res`: The oversampling to be used in frequency.
- `real_tapers`: Optional argument to decided if real or complex tapers should be provided. Default is `true` which provides real tapers.
- `tol`: Optional argument passed to the `eigs` function of `Arpack`. You likely need to play with this.

# Note about errors:
If you have an error, it is likely a convergence problem. Try setting a larger value for `tol`.
"""
function optimaltapers(
	region::Meshes.GeometryOrDomain,
	grid::CartesianGrid;
	freq_region::Meshes.GeometryOrDomain,
	ntapers::Int,
	freq_res,
	freq_downsample = nothing,
	real_tapers = true,
	tol = 0.0,
)
	@assert embeddim(grid) == embeddim(region) "The region and grid must have the same number of dimensions."
	@assert boundingbox(region) ⊆ boundingbox(grid) "The region of interest is not covered by the grid, please provide a bigger grid or smaller region."
	@assert embeddim(freq_region) == embeddim(region) "The region and freq_region must have the same number of dimensions."
	freq_res = checkfreqres(grid, freq_res)

	R = pad(downsample(pixelate_region(grid, region), freq_downsample), freq_res)
	K = pixelate_region(
		fftfreq.(freq_res, 1 ./ downsample_spacing(grid, freq_downsample)),
		freq_region,
	)
	h_oversize, λ =
		compute_eigenfunctions(R, K, ntapers; real_tapers = real_tapers, tol = tol)
	h = [
		prod(sqrt, unitless_spacing(grid)) .*
		reprocess(h_oversize[i], grid, freq_downsample) for i in eachindex(h_oversize)
	]
	return h, λ
end

function checkfreqres(grid::CartesianGrid, freq_res)
	if !(freq_res isa Int)
		@assert length(freq_res) == embeddim(grid) "freq_res must be an integer or a tuple of integers with the same length as the dimension of grid."
	end
	if all(size(grid) .≤ freq_res)
		return freq_res
	else
		@warn "freq_res is smaller than grid size, will set to the grid size but you may wish to increase this."
		return max.(size(grid), freq_res)
	end
end


"""
	reprocess(h_large, grid, freq_downsample)

Undoes the downsampling and padding.
"""
function reprocess(h_large::AbstractArray{T, D}, grid, freq_downsample) where {D, T}
	return newweight(h_large, freq_downsample) .* upsample(
		h_large[ntuple(d -> 1:downsample_size(grid, freq_downsample)[d], Val{D}())...],
		grid,
		freq_downsample,
	)
end

"""
	newweight(h, freq_downsampling)

Necessary because eigs gives normalised vector under the downsampling.
"""
newweight(::AbstractArray, ::Nothing) = 1.0
newweight(::AbstractArray{T, D}, freq_downsample::Real) where {T, D} =
	sqrt(1 / freq_downsample)^D
newweight(::AbstractArray{T, D}, freq_downsample::NTuple{D, Real}) where {T, D} =
	prod(sqrt, inv.(freq_downsample))

downsample_size(grid::CartesianGrid, ::Nothing) = size(grid)
downsample_size(grid::CartesianGrid, freq_downsample) = size(grid) .÷ freq_downsample

downsample_spacing(grid::CartesianGrid, ::Nothing) = unitless_spacing(grid)
downsample_spacing(grid::CartesianGrid, freq_downsample) =
	unitless_spacing(grid) .* freq_downsample

pixelate_region(grid::CartesianGrid, shp::Meshes.GeometryOrDomain) =
	pixelate_region(grid2side(grid), shp)
function pixelate_region(
	g::NTuple{D, AbstractArray},
	shp::Meshes.GeometryOrDomain,
) where {D}
	R = [Point(i) ∈ shp for i in Iterators.product(g...)]
	return R
end
