"""Taper concentration validation utilities."""

function taper_concentrations(taper_families, bandwidth; resolution = 100)
    concentrations = fill(
        complex(NaN),
        length(first(taper_families)),
        length(first(taper_families)),
        length(taper_families),
        length(taper_families)
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
                resolution = resolution
            )
        else
            concentrations[i, j, k, l] = conj.(concentrations[j, i, l, k])
            # because concentrations[i,j,k,l] = conj.(concentrations[j,i,l,k])
        end
    end
    return concentrations
end

function single_taper_concentration(taper1, taper2, concentration_domain; resolution)
    concentration_bb = boundingbox(concentration_domain)
    concentration_grid = CartesianGrid(
        minimum(concentration_bb),
        maximum(concentration_bb),
        dims = resolution
    )
    concentration = sum(taper_ft(taper1, centroid(concentration_grid, i)) *
                        conj(taper_ft(taper2, centroid(concentration_grid, i)))
    for i in eachindex(concentration_grid)
    if centroid(concentration_grid, i) ∈ concentration_domain)

    return prod(unitless_spacing(concentration_grid)) * concentration
end

"""
    choose_concentration_resolution(tapers_on_grids, concentration_region)

Chooses the number of wavenumbers to sample at to integrate on the concentration region.
Assumes that all tapers in each taper family are the same in some sense (the Fourier transform is sampled in the same way and interpolated).
Does not assume this across families.
Will return the lowest resolution to mitigate some effects of interpolation that can matter for fine resolutions.
"""
function choose_concentration_resolution(tapers_on_grids, concentration_region)
    bbox = boundingbox(concentration_region)
    nks = [_nk_in_box(tapers[1], bbox) for tapers in tapers_on_grids]
    concentration_res = ntuple(i -> minimum(x -> x[i], nks), Val{embeddim(bbox)}())
    if any(x -> x < 30, concentration_res)
        @warn "The resolution for computing the concentration of the tapers is very low. You may want to increase the wavenumber domain resolution used when precomputing taper Fourier transforms."
    end
    return concentration_res
end

"""
    _nk_in_box

Internal function to compute the number of wavenumbers at which the Fourier transform of a taper was evaluated in given box.
"""
function _nk_in_box(taper::ContinuousTaper, box::Box)
    return ntuple(d -> 2^63 - 1, embeddim(box))
end

function _nk_in_box(taper::Taper, box::Box)
    wavenumber = _evaluation_wavenumbers(taper)
    sides = box2sides(box)
    return ntuple(
        d -> sum(x -> sides[d][1] ≤ x ≤ sides[d][2], wavenumber[d]), embeddim(box))
end

function _evaluation_wavenumbers(taper::DiscreteTaper)
    taper.taper_ft.taper_ft.itp.ranges
end

function _evaluation_wavenumbers(taper::InterpolatedTaper)
    taper.taper_ft.taper_ft.taper_ft.itp.ranges
end
