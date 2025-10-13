# Core taper type definitions with clear discrete/continuous hierarchy.

# Abstract base types to establish the hierarchy
abstract type AbstractTaper end
abstract type DiscreteTaper <: AbstractTaper end
abstract type ContinuousTaper <: AbstractTaper end

# Helper types for discrete tapers
struct DiscreteTaperSeq{T <: AbstractArray, G <: CartesianGrid}
    taper::T
    grid::G
    function DiscreteTaperSeq(taper, grid)
        @assert size(taper)==size(grid) "The taper must have the same size as the grid."
        new{typeof(taper), typeof(grid)}(taper, grid)
    end
end

struct DiscreteTaperFT{H}
    taper_ft::H
    function DiscreteTaperFT(taper, grid, wavenumber_res)
        wavenumber_res = process_res(wavenumber_res, boundingbox(grid))
        ft_desc = fft_anydomain(
            taper, grid, wavenumber_res, 1 ./ (2 .* unitless_spacing(grid))) .*
                  prod(unitless_spacing(grid))
        wavenumber = fftshift.(fftfreq.(size(ft_desc), inv.(unitless_spacing(grid))))

        ft_interp = linear_interpolation(wavenumber, ft_desc, extrapolation_bc = Periodic())
        new{typeof(ft_interp)}(ft_interp)
    end
end

struct InterpolatedTaperFunc{H}
    taper::H
end

struct InterpolatedTaperFT{H, G <: CartesianGrid}
    taper_ft::H
    grid::G
end

# Discrete space tapers (defined on grids)
struct GridTaper{T <: AbstractArray, G <: CartesianGrid, H} <: DiscreteTaper
    values::T           # taper values on grid
    grid::G            # the grid
    fourier_transform::H   # precomputed FT for efficiency

    function GridTaper(values::T, grid::G, fourier_transform::H) where {T, G, H}
        @assert size(values)==size(grid) "Taper values must match grid size"
        new{T, G, H}(values, grid, fourier_transform)
    end
end

# Continuous space tapers derived from interpolation
struct InterpolatedTaper{I, H, G} <: ContinuousTaper
    interpolator::I         # interpolation object
    fourier_transform::H    # FT (may be discrete-based)
    source_grid::G         # original grid (for FT computation)
end

# Sin tapers
struct SinTaper{D} <: ContinuousTaper
    modes::NTuple{D, Int}    # (m1, m2, ...) mode numbers
    region::Box             # the region

    function SinTaper(modes::NTuple{D, Int}, region::Box) where {D}
        @assert all(m -> m > 0, modes) "All mode numbers must be positive"
        @assert embeddim(region)==D "Region dimension must match modes dimension"
        new{D}(modes, region)
    end
end

# Collection type
struct TaperFamily{T}
    tapers::T
end

# Collection interface
Base.eachindex(taper::TaperFamily) = eachindex(taper.tapers)
Base.length(taper::TaperFamily) = length(taper.tapers)
Base.getindex(taper::TaperFamily, i) = taper.tapers[i]

# Grid utilities
struct NoGrid end

# Utility functions for processing resolution parameters
function process_res(space_res::NTuple, bbox::Box)
    @assert embeddim(bbox)==length(space_res) "space_res must have the same dimension as the bounding box"
    space_res
end
process_res(space_res::Int, bbox::Box) = ntuple(i -> space_res, Val{embeddim(bbox)}())
process_res(space_res::Vector, bbox::Box) = process_res((Int.(space_res)...,), bbox)
process_res(space_res::Number, bbox::Box) = process_res(Int(space_res), bbox)
