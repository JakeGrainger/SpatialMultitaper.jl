# Interpolated taper construction utilities.

"""
	interp_ip_2(h::AbstractArray{T,D}, g::AbstractArray{T,D}, grid::CartesianGrid) where {T,D}

The L₂ inner product of interpolated sequences `g` and `h` recorded on grid `grid`.
"""
function L2_inner_product_interpolated(
        h::AbstractArray{T, D},
        g::AbstractArray{T, D},
        grid::CartesianGrid
) where {T, D}
    @assert size(h) == size(g) == size(grid)

    v_set = Iterators.product(ntuple(d -> 0:1, Val{D}())...)
    g_padded = centerpad(g, 1)
    h_padded = centerpad(h, 1)

    prod(unitless_spacing(grid)) * sum(
        prod((1 + (v[d] == v_prime[d])) / 6 for d in 1:D) * sum(
            gx * hx
        for (gx, hx) in zip(
            view(g_padded, ntuple(d -> v[d] .+ (1:(size(g, d) .+ 1)), Val{D}())...),
            view(h_padded, ntuple(d -> v_prime[d] .+ (1:(size(h, d) .+ 1)), Val{D}())...)
        )
        ) for v in v_set, v_prime in v_set
    )
end

function interpolate_tapers(raw_tapers, grid; wavenumber_res = size(grid))
    sides = grid2sides(grid)
    scaled_taper = raw_tapers ./
                   sqrt(L2_inner_product_interpolated(raw_tapers, raw_tapers, grid))
    taper_interp = linear_interpolation(
        sides,
        scaled_taper,
        extrapolation_bc = zero(eltype(raw_tapers))
    )
    fourier_transform = InterpolatedTaperFT(
        DiscreteTaperFT(scaled_taper, grid, wavenumber_res),
        grid
    )
    InterpolatedTaper(
        InterpolatedTaperFunc(taper_interp),
        fourier_transform,
        grid
    )
end

"""
	interpolated_taper_family(raw_tapers, grid; wavenumber_res=size(grid))

Constructs a multitaper object from multiple taper discrete impulse responses on a grid.

# Arguments
- `raw_tapers`: A `Vector` of `AbstractArray`s containing the taper impulse response.
- `grid`: The grid corresponding to the tapers.
- `wavenumber_res`: Optional named argument specifying the upsampling for computing the transfer function.

# Wavenumber sampling
- The wavenumber grid on which the transfer function will be sampled has:
	- length `wavenumber_res`.
	- interval `unitless_spacing(grid)/wavenumber_res`.
"""
function interpolated_taper_family(
        raw_tapers::Vector{<:AbstractArray},
        grid::CartesianGrid,
        concentration_region;
        wavenumber_res = size(grid)
)
    check_interpolated_tapers_cross(raw_tapers, grid; tol = 1e-2)

    interpolated_tapers = TaperFamily(
        [interpolate_tapers(taper, grid, wavenumber_res = wavenumber_res)
         for taper in raw_tapers], concentration_region
    )
    return interpolated_tapers
end

"""
    make_tapers(region; bandwidth=nothing, concentration_region=nothing, threshold=0.99, ...)
    make_tapers(region; bandwidth, threshold=0.99, ...)

A function to create tapers for a given region.

# Arguments
- `region`: The region to create the tapers on.
- `bandwidth`: The bandwidth of the tapers (creates a Ball concentration region)
- `concentration_region`: Direct specification of the concentration region in wavenumber space
- `threshold`: The threshold for the number of tapers to use.
- `space_res`: The resolution of the grid to create the tapers on.
- `wavenumber_res`: The resolution of the wavenumber grid to compute the transfer function on.
- `wavenumber_generate_res`: The resolution of the wavenumber grid to generate the tapers on.
- `ntapers_max`: The maximum number of tapers to generate.

Either `bandwidth` or `concentration_region` must be provided, but not both.

As a rough rule of thumb, setting the bandwidth to an integer multiple (say `x`) of a side length of the bounding square of a region will
produce `x^d` tapers (note this is a very crude approximation, but you will find out how many good tapers you can make when you create them).
"""
function make_tapers(
        region;
        bandwidth = nothing,
        concentration_region = nothing,
        threshold = 0.99,
        space_res = 100,
        wavenumber_res = 500,
        wavenumber_generate_res = 500,
        ntapers_max = 30
)
    # Validate arguments
    if bandwidth !== nothing && concentration_region !== nothing
        throw(ArgumentError("Cannot specify both bandwidth and concentration_region"))
    elseif bandwidth === nothing && concentration_region === nothing
        throw(ArgumentError("Must specify either bandwidth or concentration_region"))
    end

    # Create concentration region from bandwidth if needed
    if bandwidth !== nothing
        D = embeddim(region)
        concentration_region = Ball(Point(ntuple(_ -> 0.0, D)), bandwidth)
    end

    bbox = boundingbox(region)
    space_res_processed = process_res(space_res, bbox)
    wavenumber_res_processed = process_res(wavenumber_res, bbox)
    wavenumber_generate_res_processed = process_res(wavenumber_generate_res, bbox)

    grid = CartesianGrid(bbox.min, bbox.max, dims = space_res_processed)
    # TODO: need to do some rounding here to make this more reliable
    @info "Computing tapers..."
    h, λ = optimaltapers(
        region,
        grid,
        wavenumber_region = concentration_region,
        ntapers = Int(ntapers_max),
        wavenumber_res = wavenumber_generate_res_processed,
        wavenumber_downsample = nothing,
        tol = 1e-8,
        check_grid = false
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
    tapers = interpolated_taper_family(
        h[1:ntapers], grid, concentration_region, wavenumber_res = wavenumber_res_processed)
    @info "...made Fourier transform."

    return tapers
end
