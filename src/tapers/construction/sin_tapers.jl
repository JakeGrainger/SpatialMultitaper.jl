# Sin taper family construction utilities.

"""
    sin_taper_family(ntapers, region::Box)

Construct the sin tapers for a given region.

# Arguments:
- `ntapers`: The number of tapers per dimension (Int or Tuple of Ints).
- `region`: The region to construct the tapers on. Must be a box.

# Background:
The sin tapers on a unit interval are defined as:
```math
h_m(x) = \\sqrt{2} * \\sin(πxm)
```
"""
function sin_taper_family(ntapers::NTuple{D, Int}, region::Box) where {D}
    @assert embeddim(region)==D "Region dimension must match ntapers dimension"
    tapers = [SinTaper(modes, region)
              for modes in Iterators.product(ntuple(d -> 1:ntapers[d], D)...)]
    TaperFamily(vec(tapers))
end

# Convenience method for uniform ntapers across all dimensions
function sin_taper_family(ntapers::Int, region::Box{D}) where {D}
    sin_taper_family(ntuple(_ -> ntapers, D), region)
end

# Legacy method name for backward compatibility
make_sin_taper(modes, region) = SinTaper(modes, region)

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
