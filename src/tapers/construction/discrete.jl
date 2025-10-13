# Discrete taper construction utilities.

function discretetaper(
        taper::AbstractArray,
        grid::CartesianGrid;
        wavenumber_res = size(grid),
        rescale = true
)
    @assert size(taper) == size(grid)
    taper_scaled = if rescale
        taper .* sqrt(prod(unitless_spacing(grid)) / sum(abs2, taper))
    else
        taper
    end

    disc_taper = DiscreteTaperSeq(taper_scaled, grid)
    ft_interp = DiscreteTaperFT(taper_scaled, grid, wavenumber_res)

    Taper(disc_taper, ft_interp)
end
