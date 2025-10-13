"""Taper normalization validation utilities."""

function tapers_normalisations(taper_families)
    [[single_taper_normalisations(taper_families[j][i])
      for i in eachindex(taper_families[j])] for j in eachindex(taper_families)]
end

function single_taper_normalisations(::ContinuousTaper)
    return 1.0 # assumed normalised in this case
end

function single_taper_normalisations(taper::DiscreteTaper)
    return prod(unitless_spacing(taper.taper.grid)) * sum(abs2, taper.taper.taper)
end
