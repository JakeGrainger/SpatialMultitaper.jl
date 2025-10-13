# Core taper type definitions.

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

struct Taper{F, H}
    taper::F
    taper_ft::H
end

# Type aliases for clarity
const DiscreteTaper = Taper{<:DiscreteTaperSeq, <:DiscreteTaperFT}
const InterpolatedTaper = Taper{<:InterpolatedTaperFunc, <:InterpolatedTaperFT}
const ContinuousTaper = Union{InterpolatedTaper, Taper{<:Function, <:Function}}

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
