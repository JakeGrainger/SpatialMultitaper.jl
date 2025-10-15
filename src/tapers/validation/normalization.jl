# Taper normalization validation utilities.

function tapers_normalisations(taper_families)
    [[single_taper_normalisations(taper_families[j][i])
      for i in eachindex(taper_families[j])] for j in eachindex(taper_families)]
end

function single_taper_normalisations(::ContinuousTaper)
    return 1.0 # assumed normalised in this case
end

function single_taper_normalisations(taper::GridTaper)
    return prod(unitless_spacing(taper.grid)) * sum(abs2, taper.values)
end

# Legacy methods for backward compatibility
single_taper_normalisations(taper::DiscreteTaper) = single_taper_normalisations(taper)
