## types

struct DiscreteTaperSeq{T<:AbstractArray,G<:CartesianGrid}
    taper::T
    grid::G
    function DiscreteTaperSeq(taper, grid)
        @assert size(taper) == size(grid) "The taper must have the same size as the grid."
        new{typeof(taper),typeof(grid)}(taper, grid)
    end
end

struct DiscreteTaperTF{H}
    taper_ft::H
    function DiscreteTaperTF(taper, grid, freq_res)
        ft_desc = fftshift(fft(pad(taper, freq_res))) .* prod(unitless_spacing(grid))
        freq = fftshift.(fftfreq.(size(ft_desc), inv.(unitless_spacing(grid))))
        ft_interp = linear_interpolation(freq, ft_desc, extrapolation_bc = Periodic())
        new{typeof(ft_interp)}(ft_interp)
    end
end

struct InterpolatedTaperFunc{H}
    taper::H
end

struct InterpolatedTaperFT{H,G<:CartesianGrid}
    taper_ft::H
    grid::G
end

struct Taper{F,H}
    taper::F
    taper_ft::H
end
const DiscreteTaper = Taper{<:DiscreteTaperSeq,<:DiscreteTaperTF}
const InterpolatedTaper = Taper{<:InterpolatedTaperFunc,<:InterpolatedTaperFT}
const ContinuousTaper = Union{InterpolatedTaper,Taper{<:Function,<:Function}}

struct TaperFamily{M,T<:NTuple{M,Taper}}
    tapers::T
end
Base.eachindex(taper::TaperFamily) = eachindex(taper.tapers)
Base.length(taper::TaperFamily) = length(taper.tapers)
Base.getindex(taper::TaperFamily, i) = taper.tapers[i]

## functions

function discretetaper(taper::AbstractArray, grid::CartesianGrid; freq_res = size(grid))
    @assert size(taper) == size(grid)
    taper_scaled = taper ./ sqrt(sum(abs2, taper) * prod(unitless_spacing(grid)))

    disc_taper = DiscreteTaperSeq(taper_scaled, grid)
    ft_interp = DiscreteTaperTF(taper_scaled, grid, freq_res)

    Taper(disc_taper, ft_interp)
end

"""
	interp_ip_2(h::AbstractArray{T,D}, g::AbstractArray{T,D}, grid::CartesianGrid) where {T,D}

The L₂ inner product of interpolated sequences `g` and `h` recorded on grid `grid`.
"""
function L2_inner_product_interpolated(
    h::AbstractArray{T,D},
    g::AbstractArray{T,D},
    grid::CartesianGrid,
) where {T,D}
    @assert size(h) == size(g) == size(grid)

    v_set = Iterators.product(ntuple(d -> 0:1, Val{D}())...)
    g_padded = pad(g, size(g) .+ 1)
    h_padded = pad(h, size(h) .+ 1)

    prod(unitless_spacing(grid)) * sum(
        prod((-1)^(v[d] + v_prime[d]) / 3 + v[d] * v_prime[d] for d = 1:D) * sum(
            hx * gx for (hx, gx) in zip(
                view(g_padded, ntuple(d -> v[d]+1:size(g, d), Val{D}())...),
                view(h_padded, ntuple(d -> v[d]+1:size(g, d), Val{D}())...),
            )
        ) for v in v_set, v_prime in v_set
    )
end

function Interpolations.interpolate(taper::DiscreteTaper)
    sides = grid2side(taper.taper.grid)
    scaled_taper =
        taper.taper.taper ./ sqrt(
            L2_inner_product_interpolated(
                taper.taper.taper,
                taper.taper.taper,
                taper.taper.grid,
            ),
        )
    taper_interp = linear_interpolation(
        sides,
        scaled_taper,
        extrapolation_bc = zero(eltype(taper.taper.taper)),
    )
    Taper(
        InterpolatedTaperFunc(taper_interp),
        InterpolatedTaperFT(taper.taper_ft, taper.taper.grid),
    )
end

## taper evaluation
(taper::Taper)(x...) = taper(x)
(taper::Taper)(x::Point) = taper(unitless_coords(x))
function (taper::Taper)(x)
    taper.taper(x)
end
taper_ft(taper::Taper, x...) = taper_ft(taper, x)
taper_ft(taper::Taper, x::Point) = taper_ft(taper, unitless_coords(x))
function taper_ft(taper::Taper, x)
    taper.taper_ft(x)
end

(f::InterpolatedTaperFunc)(x::NTuple{D,Real}) where {D,Real} = f.taper(x...)
(f::InterpolatedTaperFunc)(x::NTuple{1,Real}) = f.taper(x[1])
(f::InterpolatedTaperFunc)(x::NTuple{2,Real}) = f.taper(x[1], x[2])
(f::InterpolatedTaperFunc)(x::NTuple{3,Real}) = f.taper(x[1], x[2], x[3])

(f::InterpolatedTaperFT)(k::NTuple{D,Real}) where {D,Real} =
    f.taper_ft(k) * prod(sinc(π * unitless_spacing(f.grid)[d] * k[d])^2 for d = 1:D)

(f::DiscreteTaperTF)(x::NTuple{D,Real}) where {D,Real} = f.taper_ft(x...)
(f::DiscreteTaperTF)(x::NTuple{1,Real}) = f.taper_ft(x[1])
(f::DiscreteTaperTF)(x::NTuple{2,Real}) = f.taper_ft(x[1], x[2])
(f::DiscreteTaperTF)(x::NTuple{3,Real}) = f.taper_ft(x[1], x[2], x[3])

## show methods
Base.show(io::IO, ::MIME"text/plain", taper::ContinuousTaper) =
    print(io, "A continuous taper function.")
Base.show(io::IO, ::MIME"text/plain", taper::DiscreteTaper) =
    print(io, "A discrete taper function.")
Base.show(io::IO, ::MIME"text/plain", taper::InterpolatedTaper) =
    print(io, "A continuous interpolated taper function.")
Base.show(io::IO, ::MIME"text/plain", taper::TaperFamily) =
    print(io, "A family of $(length(taper)) taper functions.")

##

"""
	interpolated_taper_family(raw_tapers, grid; freq_res=size(grid))

Constructs a multitaper object from multiple taper discrete impulse responses on a grid. 

# Arguments
- `raw_tapers`: A `Vector` of `AbstractArray`s containing the taper impulse response.
- `grid`: The grid corresponding to the tapers.
- `freq_res`: Optional named arguemnt specifying the upsampling for computing the transfer function.

# Frequency sampling
- The frequency grid on which the transfer function will be sampled has:
	- length `freq_res`.
	- interval `unitless_spacing(grid)/freq_res`.
"""
function interpolated_taper_family(
    raw_tapers::Vector{<:AbstractArray},
    grid::CartesianGrid;
    freq_res = size(grid),
)
    interpolated_tapers = TaperFamily(
        ntuple(
            d -> interpolate(discretetaper(raw_tapers[d], grid, freq_res = freq_res)),
            length(raw_tapers),
        ),
    )
    # checks orthogonality of interpolated tapers
    max_cross = maximum(
        abs(L2_inner_product_interpolated(raw_tapers[i], raw_tapers[j], grid)) for
        i = 1:length(interpolated_tapers) for j = i+1:length(interpolated_tapers)
    )
    if max_cross > 1e-2
        max_discrete_cross = maximum(
            abs(sum(h * g for (h, g) in zip(raw_tapers[i], raw_tapers[j]))) for
            i = 1:length(interpolated_tapers) for j = i+1:length(interpolated_tapers)
        )
        @warn "The tapers are not orthogonal. The maximum cross-correlation is $max_cross. You may wish to increase the resolution of the grid used for the discrete tapers, or check that the discrete tapers are orthogonal (the discrete tapers has a max-cross L2 norm of $max_discrete_cross)."
    end

    return interpolated_tapers
end


"""
	make_sin_tapers(lbw, region::Box, gridsize)

Construct the sin tapers for a given region and gridsize.

# Arguments:
- `ntapers`: The number of tapers.
- `region`: The region to construct the tapers on. Must be a box.

# Background:
The sin tapers on a unit interval are defined as:
```math
h_m(x) = \\sqrt{2} * \\sin(πxm)
```
"""
function sin_taper_family(ntapers, region::Box)
    tapers =
        [make_sin_taper(d, region) for d in Iterators.ProductIterator(range.(1, ntapers))]
    TaperFamily(ntuple(d -> tapers[d], length(tapers)))
end

function make_sin_taper(m, region)
    region_sides = getfield.(sides(region), :val)
    region_start = unitless_coords(minimum(region))
    phase = exp(-2π * 1im * sum(region_start))
    function _single_sin_taper(x)
        prod(
            sin_taper(x[d] - region_start[d], m[d], region_sides[d]) for
            d = 1:embeddim(region)
        )
    end
    function _single_sin_taper_ft(k)
        prod(sin_ft(k[d], m[d], region_sides[d]) for d = 1:embeddim(region)) * phase
    end
    return Taper(_single_sin_taper, _single_sin_taper_ft)
end

sin_taper_base(x, m) = sqrt(2) * sin(π * x * m) * (0 ≤ x ≤ 1)
sin_taper(x, m, l) = sin_taper_base(x / l, m) / sqrt(l)

function sin_ft_base(k::Real, m::Int; tol = 1e-10)
    if abs(2k - m) < tol
        return -1im / sqrt(2)
    elseif abs(2k + m) < tol
        return 1im / sqrt(2)
    else
        return exp(-1im * π * (k - (m - 1) / 2)) / sqrt(2) * m * sin(π * k - π * m / 2) /
               (π * (k^2 - (m / 2)^2))
    end
end
sin_ft(k, m, l) = sin_ft_base(k * l, m) * sqrt(l)

"""
    make_tapers(region; bandwidth, threshold = 0.99, space_res = 100, freq_res = 500, freq_generate_res = 500, ntapers_max = 30)

A function to create tapers for a given region.
	
# Arguments
- `region`: The region to create the tapers on.
- `bandwidth`: The bandwidth of the tapers
- `threshold`: The threshold for the number of tapers to use.
- `space_res`: The resolution of the grid to create the tapers on.
- `freq_res`: The resolution of the frequency grid to compute the transfer function on.
- `freq_generate_res`: The resolution of the frequency grid to generate the tapers on.
- `ntapers_max`: The maximum number of tapers to generate.

As a rough rule of thumb, setting the bandwidth to an integer multiple (say `x`) of a side length of the bounding square of a region will 
produce `x^d` tapers (note this is a very crude approximation, but you will find out how many good tapers you can make when you create them).
"""
function make_tapers(
    region;
    bandwidth,
    threshold = 0.99,
    space_res = 100,
    freq_res = 500,
    freq_generate_res = 500,
    ntapers_max = 30,
)
    bbox = boundingbox(region)
    space_res_processed = process_res(space_res, bbox)
    freq_res_processed = process_res(freq_res, bbox)
    freq_generate_res_processed = process_res(freq_generate_res, bbox)

    grid = CartesianGrid(bbox.min, bbox.max, dims = space_res_processed)

    freq_region = Ball(Meshes.Point(0, 0), bandwidth)
    @info "Computing tapers..."
    h, λ = optimaltapers(
        region,
        grid,
        freq_region = freq_region,
        ntapers = Int(ntapers_max),
        freq_res = freq_generate_res_processed,
        freq_downsample = nothing,
        tol = 1e-8,
    )
    @info "...tapers computed."
    ## process tapers
    @info "Finding number of tapers to use..."
    ntaper_index = findfirst(real.(λ) .< threshold)
    ntapers = isnothing(ntaper_index) ? length(λ) : ntaper_index - 1
    @assert ntapers ∈ eachindex(λ)
    if ntapers == length(λ)
        @warn "All tapers are used, consider increasing the threshold or the number of pregenerated tapers"
    end
    @info "...number of tapers to use: $ntapers out of $(length(λ)) available."

    @info "Making transfer function..."
    tapers = interpolated_taper_family(h[1:ntapers], grid, freq_res = freq_res_processed)
    @info "...made transfer function."

    return tapers
end

function process_res(space_res::NTuple, bbox::Box)
    @assert embeddim(bbox) == length(space_res) "space_res must have the same dimension as the bounding box"
    space_res
end
process_res(space_res::Int, bbox::Box) = ntuple(i -> space_res, Val{embeddim(bbox)}())
process_res(space_res::Vector, bbox::Box) = process_res((Int.(space_res)...,), bbox)
process_res(space_res::Number, bbox::Box) = process_res(Int(space_res), bbox)