# Discrete taper construction utilities.

function gridtaper(
        values::AbstractArray,
        grid::CartesianGrid;
        wavenumber_res = size(grid),
        rescale = true
)
    @assert size(values) == size(grid)
    values_scaled = if rescale
        values .* sqrt(prod(unitless_spacing(grid)) / sum(abs2, values))
    else
        values
    end

    fourier_transform = DiscreteTaperFT(values_scaled, grid, wavenumber_res)
    GridTaper(values_scaled, grid, fourier_transform)
end

# Legacy function name for backward compatibility
function discretetaper(
        taper::AbstractArray,
        grid::CartesianGrid;
        wavenumber_res = size(grid),
        rescale = true
)
    gridtaper(taper, grid; wavenumber_res = wavenumber_res, rescale = rescale)
end
