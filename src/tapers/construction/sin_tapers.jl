# Sin taper family construction utilities.

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
    tapers = [make_sin_taper(d, region)
              for d in Iterators.ProductIterator(range.(1, ntapers))]
    TaperFamily(tapers[:])
end

function make_sin_taper(m, region)
    region_sides = getfield.(sides(region), :val)
    region_start = unitless_coords(minimum(region))
    phase = exp(-2π * 1im * sum(region_start))
    function _single_sin_taper(x)
        prod(
            sin_taper(x[d] - region_start[d], m[d], region_sides[d])
        for
        d in 1:embeddim(region)
        )
    end
    function _single_sin_taper_ft(k)
        prod(sin_ft(k[d], m[d], region_sides[d]) for d in 1:embeddim(region)) * phase
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
