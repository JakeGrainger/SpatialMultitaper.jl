## types

struct DiscreteTaperSeq{T<:AbstractArray,G<:CartesianGrid}
    taper::T
    grid::G
    function DiscreteTaperSeq(taper, grid)
        @assert size(taper) == size(grid) "The taper must have the same size as the grid."
        new{typeof(taper),typeof(grid)}(taper, grid)
    end
end

struct DiscreteTaperFT{H}
    taper_ft::H
    function DiscreteTaperFT(taper, grid, freq_res)
        ft_desc = fftshift(fft(pad(taper, freq_res))) .* prod(unitless_spacing(grid))
        freq = fftshift.(fftfreq.(size(ft_desc), inv.(unitless_spacing(grid))))
        ft_interp = cubic_spline_interpolation(freq, ft_desc, extrapolation_bc = Periodic())
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
const DiscreteTaper = Taper{<:DiscreteTaperSeq,<:DiscreteTaperFT}
const InterpolatedTaper = Taper{<:InterpolatedTaperFunc,<:InterpolatedTaperFT}
const ContinuousTaper = Union{InterpolatedTaper,Taper{<:Function,<:Function}}

struct TaperFamily{M,T<:NTuple{M,Taper}}
    tapers::T
end
Base.eachindex(taper::TaperFamily) = eachindex(taper.tapers)
Base.length(taper::TaperFamily) = length(taper.tapers)
Base.getindex(taper::TaperFamily, i) = taper.tapers[i]

## functions
function discretetaper(
    taper::AbstractArray,
    grid::CartesianGrid;
    freq_res = size(grid),
    rescale = true,
)
    @assert size(taper) == size(grid)
    taper_scaled = if rescale
        taper .* sqrt(prod(unitless_spacing(grid)) / sum(abs2, taper))
    else
        taper
    end

    disc_taper = DiscreteTaperSeq(taper_scaled, grid)
    ft_interp = DiscreteTaperFT(taper_scaled, grid, freq_res)

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
    g_padded = centerpad(g, 1)
    h_padded = centerpad(h, 1)

    prod(unitless_spacing(grid)) * sum(
        prod((1 + (v[d] == v_prime[d])) / 6 for d = 1:D) * sum(
            gx * hx for (gx, hx) in zip(
                view(g_padded, ntuple(d -> v[d] .+ (1:size(g, d).+1), Val{D}())...),
                view(h_padded, ntuple(d -> v_prime[d] .+ (1:size(h, d).+1), Val{D}())...),
            )
        ) for v in v_set, v_prime in v_set
    )
end

function Interpolations.interpolate(taper::DiscreteTaper)
    Interpolations.interpolate
end

function interpolate_tapers(raw_tapers, grid; freq_res = size(grid))
    sides = grid2side(grid)
    scaled_taper =
        raw_tapers ./ sqrt(L2_inner_product_interpolated(raw_tapers, raw_tapers, grid))
    taper_interp = linear_interpolation(
        sides,
        scaled_taper,
        extrapolation_bc = zero(eltype(raw_tapers)),
    )
    Taper(
        InterpolatedTaperFunc(taper_interp),
        InterpolatedTaperFT(DiscreteTaperFT(scaled_taper, grid, freq_res), grid),
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
    f.taper_ft(k) * prod(sinc(unitless_spacing(f.grid)[d] * k[d])^2 for d = 1:D)

(f::DiscreteTaperFT)(x::NTuple{D,Real}) where {D,Real} = f.taper_ft(x...)
(f::DiscreteTaperFT)(x::NTuple{1,Real}) = f.taper_ft(x[1])
(f::DiscreteTaperFT)(x::NTuple{2,Real}) = f.taper_ft(x[1], x[2])
(f::DiscreteTaperFT)(x::NTuple{3,Real}) = f.taper_ft(x[1], x[2], x[3])

## show methods
Base.show(io::IO, ::MIME"text/plain", taper::ContinuousTaper) =
    print(io, "A continuous taper function.")
Base.show(io::IO, ::MIME"text/plain", taper::DiscreteTaper) =
    print(io, "A discrete taper function.")
Base.show(io::IO, ::MIME"text/plain", taper::InterpolatedTaper) =
    print(io, "A continuous interpolated taper function.")
Base.show(io::IO, ::MIME"text/plain", taper::TaperFamily) =
    print(io, "A family of $(length(taper)) taper functions.")

struct NoGrid end
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
function tapers_on_grid(tapers::TaperFamily, ::NoGrid; freq_res = 500) # currently ignore freq_res
    return tapers
end
function tapers_on_grid(
    tapers::TaperFamily{M},
    grid::CartesianGrid;
    freq_res = 500,
) where {M}
    return TaperFamily(
        ntuple(i -> single_taper_on_grid(tapers[i], grid; freq_res = freq_res), Val{M}()),
    )
end
function single_taper_on_grid(taper::ContinuousTaper, grid::CartesianGrid; freq_res = 500)
    taper_evaluated = [taper(x) for x in Iterators.ProductIterator(grid2side(grid))]
    return discretetaper(taper_evaluated, grid, freq_res = freq_res, rescale = false) # don't rescale as taper transform gets rescaled later and we account for this in the norm
end

function tapers_normalisations(taper_families)
    [
        [
            single_taper_normalisations(taper_families[j][i]) for
            i in eachindex(taper_families[j])
        ] for j in eachindex(taper_families)
    ]
end
function single_taper_normalisations(::ContinuousTaper)
    1.0 # assumed normalised in this case
end
function single_taper_normalisations(taper::DiscreteTaper)
    prod(unitless_spacing(taper.taper.grid)) * sum(abs2, taper.taper.taper)
end

function taper_concentrations(taper_families, bandwidth; resolution = 100)
    concentrations = fill(
        complex(NaN),
        length(first(taper_families)),
        length(first(taper_families)),
        length(taper_families),
        length(taper_families),
    )
    for l in axes(concentrations, 4),
        k in axes(concentrations, 3),
        j in axes(concentrations, 2),
        i in axes(concentrations, 1)

        if isnan(concentrations[j, i, l, k])
            concentrations[i, j, k, l] = single_taper_concentration(
                taper_families[k][i],
                taper_families[l][j],
                bandwidth,
                resolution = resolution,
            )
        else
            concentrations[i, j, k, l] = conj.(concentrations[j, i, l, k])
            # because concentrations[i,j,k,l] = conj.(concentrations[j,i,l,k])
        end
    end
    return concentrations
end
function single_taper_concentration(taper1, taper2, concentration_domain; resolution = 100)
    concentration_bb = boundingbox(concentration_domain)
    concentration_grid = CartesianGrid(
        minimum(concentration_bb),
        maximum(concentration_bb),
        dims = ntuple(d -> resolution, Val{embeddim(concentration_bb)}()),
    )
    concentration =
        prod(unitless_spacing(concentration_grid)) * sum(
            taper_ft(taper1, x) * conj(taper_ft(taper2, x)) for
            x in Iterators.ProductIterator(grid2side(concentration_grid)) if
            Point(x) ∈ concentration_domain
        )
    return concentration
end


"""
    check_tapers_for_data(data, tapers, bandwidth)

This function checks if the tapers are suitable for the provided data.
In particular, it checks that the tapers work for the type of grids provided in the data.
This assumes that the base tapers have an L2 norm of 1 and are continuous.
We use the term `grid` here to refer to the sampling mechanism of a process.
So it also refer to continuously recorded data (which is the `NoGrid` type).

If you have a problem, you can use `taper_checks` to check the tapers directly.
"""
function check_tapers_for_data(data, tapers, bandwidth; freq_res = 500, tol = 1e-2)
    normalisations, concentrations =
        taper_checks(data, tapers, bandwidth, freq_res = freq_res)
    concentrations_diag = [
        concentrations[i, i, k, l] for i in axes(concentrations, 1),
        k in axes(concentrations, 3), l in axes(concentrations, 4)
    ]
    concentrations_off_diag = [
        i == j ? zero(eltype(concentrations)) : concentrations[i, j, k, l] for
        i in axes(concentrations, 1), j in axes(concentrations, 2),
        k in axes(concentrations, 3), l in axes(concentrations, 4)
    ]
    @assert all(x -> abs(x - 1) < tol, normalisations) "All tapers must have an L2 norm of 1, but found largest difference of $(maximum(x->abs(x-1), normalisations))"
    @assert all(x -> abs(x - 1) < tol, concentrations_diag) "All tapers must have a concentration of 1, but found largest difference of $(maximum(x->abs(x-1), concentrations_diag))"
    @assert all(x -> abs(x) < tol, concentrations_off_diag) "All tapers must have no cross-concentration, but found largest difference of $(maximum(x->abs(x), concentrations_off_diag))"
end
function taper_checks(data, tapers, bandwidth; freq_res = 500, conc_res = 100)
    @assert bandwidth > 0 "Bandwidth must be positive"

    observational_types = unique(get_grid.(data))
    tapers_on_grids = tapers_on_grid.(Ref(tapers), observational_types; freq_res = freq_res)
    concentration_region = Ball(Point(ntuple(i -> 0.0, _getdims(data[1]))), bandwidth)
    normalisations = tapers_normalisations(tapers_on_grids) # checks on each family
    concentrations =
        taper_concentrations(tapers_on_grids, concentration_region; resolution = conc_res) # has to check across the sampled and non sampled versions
    return normalisations, concentrations
end

## specific taper families
"""
	interpolated_taper_family(raw_tapers, grid; freq_res=size(grid))

Constructs a multitaper object from multiple taper discrete impulse responses on a grid. 

# Arguments
- `raw_tapers`: A `Vector` of `AbstractArray`s containing the taper impulse response.
- `grid`: The grid corresponding to the tapers.
- `freq_res`: Optional named argument specifying the upsampling for computing the transfer function.

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
    check_interpolated_tapers_cross(raw_tapers, grid; tol = 1e-2)

    interpolated_tapers = TaperFamily(
        ntuple(
            d -> interpolate_tapers(raw_tapers[d], grid, freq_res = freq_res),
            length(raw_tapers),
        ),
    )
    return interpolated_tapers
end

"""
    check_interpolated_tapers_cross(raw_tapers, grid; tol = 1e-2)

Checks how much the interpolation changes the orthogonality of the tapers. If the maximum cross-correlation is greater than `1e-2`, it will warn the user.
"""
function check_interpolated_tapers_cross(raw_tapers, grid; tol = 1e-2)
    max_cross = maximum(
        abs(L2_inner_product_interpolated(raw_tapers[i], raw_tapers[j], grid)) for
        i in eachindex(raw_tapers) for j = i+1:length(raw_tapers)
    )
    if max_cross > tol
        max_discrete_cross = maximum(
            abs(sum(h * g for (h, g) in zip(raw_tapers[i], raw_tapers[j]))) for
            i in eachindex(raw_tapers) for j = i+1:length(raw_tapers)
        )
        @warn "The tapers are significantly not orthogonal. The maximum cross-correlation is $max_cross. You may wish to increase the resolution of the grid used for the discrete tapers, or check that the discrete tapers are orthogonal (the discrete tapers has a max-cross L2 norm of $max_discrete_cross)."
    end
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

    @info "Making Fourier transform of tapers..."
    tapers = interpolated_taper_family(h[1:ntapers], grid, freq_res = freq_res_processed)
    @info "...made Fourier transform."

    return tapers, λ[1:ntapers]
end

function process_res(space_res::NTuple, bbox::Box)
    @assert embeddim(bbox) == length(space_res) "space_res must have the same dimension as the bounding box"
    space_res
end
process_res(space_res::Int, bbox::Box) = ntuple(i -> space_res, Val{embeddim(bbox)}())
process_res(space_res::Vector, bbox::Box) = process_res((Int.(space_res)...,), bbox)
process_res(space_res::Number, bbox::Box) = process_res(Int(space_res), bbox)