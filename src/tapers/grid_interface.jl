# Grid conversion and interface utilities.

"""
    get_grid(data)

Gets the grid associated with `data`. If no grid is associated, it returns `NoGrid()`, i.e. if `data is PointSet` or `domain(data) isa PointSet`.
"""
get_grid(data::GeoTable) = get_grid(domain(data))
get_grid(::PointSet) = NoGrid()
get_grid(grid::CartesianGrid) = grid

"""
    tapers_on_grid(tapers::TaperFamily, grid)

Converts continuous tapers to a discrete taper family on the provided grid.
If the grid is `NoGrid`, it returns the tapers unchanged.
This is used for checking if the tapers are suitable for the combination of grids present in the data.
"""
function tapers_on_grid(tapers::TaperFamily, ::NoGrid; wavenumber_res = 500) # currently ignore wavenumber_res
    return tapers
end

function tapers_on_grid(tapers::TaperFamily, grid::CartesianGrid; wavenumber_res = 500)
    return TaperFamily([single_taper_on_grid(
                            tapers[i], grid; wavenumber_res = wavenumber_res)
                        for i in eachindex(tapers)])
end

function single_taper_on_grid(
        taper::ContinuousTaper, grid::CartesianGrid; wavenumber_res = 500)
    taper_evaluated = reshape(
        map(i -> taper(centroid(grid, i)), eachindex(grid)), size(grid))
    return discretetaper(
        taper_evaluated, grid, wavenumber_res = wavenumber_res, rescale = false) # don't rescale as taper transform gets rescaled later and we account for this in the norm
end
