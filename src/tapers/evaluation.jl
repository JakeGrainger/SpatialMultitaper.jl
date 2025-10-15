# Taper evaluation and function call methods.

# Main taper evaluation interface - works for all AbstractTaper subtypes
(taper::AbstractTaper)(x...) = taper(x)
(taper::AbstractTaper)(x::Point) = taper(unitless_coords(x))

# GridTaper evaluation
function (taper::GridTaper)(x::NTuple{D, Real}) where {D}
    # For discrete tapers, we would need interpolation to evaluate at arbitrary points
    # For now, throw an error suggesting to use InterpolatedTaper instead
    throw(ArgumentError("GridTaper cannot be evaluated at arbitrary points. Convert to InterpolatedTaper first."))
end

# InterpolatedTaper evaluation
function (taper::InterpolatedTaper)(x::NTuple{D, Real}) where {D}
    taper.interpolator(x)
end

# SinTaper evaluation
function (taper::SinTaper{D})(x::NTuple{D, Real}) where {D}
    region_sides = getfield.(sides(taper.region), :val)
    region_start = unitless_coords(minimum(taper.region))

    prod(
        sin_taper(x[d] - region_start[d], taper.modes[d], region_sides[d])
    for d in 1:D
    )
end

# Fourier transform evaluation interface
taper_ft(taper::AbstractTaper, x...) = taper_ft(taper, x)
taper_ft(taper::AbstractTaper, x::Point) = taper_ft(taper, unitless_coords(x))

# GridTaper FT evaluation
function taper_ft(taper::GridTaper, k::NTuple{D, Real}) where {D}
    taper.fourier_transform(k)
end

# InterpolatedTaper FT evaluation
function taper_ft(taper::InterpolatedTaper, k::NTuple{D, Real}) where {D}
    taper.fourier_transform(k) *
    prod(sinc(unitless_spacing(taper.source_grid)[d] * k[d])^2 for d in 1:D)
end

# SinTaper FT evaluation - analytical formula
function taper_ft(taper::SinTaper{D}, k::NTuple{D, Real}) where {D}
    region_sides = getfield.(sides(taper.region), :val)
    region_start = unitless_coords(minimum(taper.region))
    phase = exp(-2π * 1im * sum(region_start))

    prod(sin_ft(k[d], taper.modes[d], region_sides[d]) for d in 1:D) * phase
end

# Interpolated taper function evaluation
(f::InterpolatedTaperFunc)(x::NTuple{D, Real}) where {D} = f.taper(x...)
(f::InterpolatedTaperFunc)(x::NTuple{1, Real}) = f.taper(x[1])
(f::InterpolatedTaperFunc)(x::NTuple{2, Real}) = f.taper(x[1], x[2])
(f::InterpolatedTaperFunc)(x::NTuple{3, Real}) = f.taper(x[1], x[2], x[3])

# Interpolated taper FT evaluation
function (f::InterpolatedTaperFT)(k::NTuple{D, Real}) where {D}
    f.taper_ft(k) * prod(sinc(unitless_spacing(f.grid)[d] * k[d])^2 for d in 1:D)
end

# Discrete taper FT evaluation
(f::DiscreteTaperFT)(x::NTuple{D, Real}) where {D} = f.taper_ft(x...)
(f::DiscreteTaperFT)(x::NTuple{1, Real}) = f.taper_ft(x[1])
(f::DiscreteTaperFT)(x::NTuple{2, Real}) = f.taper_ft(x[1], x[2])
(f::DiscreteTaperFT)(x::NTuple{3, Real}) = f.taper_ft(x[1], x[2], x[3])
