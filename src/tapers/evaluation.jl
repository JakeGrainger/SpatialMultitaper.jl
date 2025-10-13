# Taper evaluation and function call methods.

# Main taper evaluation interface
(taper::Taper)(x...) = taper(x)
(taper::Taper)(x::Point) = taper(unitless_coords(x))
function (taper::Taper)(x)
    taper.taper(x)
end

# Fourier transform evaluation
taper_ft(taper::Taper, x...) = taper_ft(taper, x)
taper_ft(taper::Taper, x::Point) = taper_ft(taper, unitless_coords(x))
function taper_ft(taper::Taper, x)
    taper.taper_ft(x)
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
