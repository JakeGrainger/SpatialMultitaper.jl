"""Taper orthogonality validation utilities."""

"""
    check_interpolated_tapers_cross(raw_tapers, grid; tol = 1e-2)

Checks how much the interpolation changes the orthogonality of the tapers. If the maximum cross-correlation is greater than `1e-2`, it will warn the user.
"""
function check_interpolated_tapers_cross(raw_tapers, grid; tol = 1e-2)
    max_cross = maximum(
        abs(L2_inner_product_interpolated(raw_tapers[i], raw_tapers[j], grid))
    for
    i in eachindex(raw_tapers) for j in (i + 1):length(raw_tapers)
    )
    if max_cross > tol
        max_discrete_cross = maximum(
            abs(sum(h * g for (h, g) in zip(raw_tapers[i], raw_tapers[j])))
        for
        i in eachindex(raw_tapers) for j in (i + 1):length(raw_tapers)
        )
        @warn "The tapers are significantly not orthogonal. The maximum cross-correlation is $max_cross. You may wish to increase the resolution of the grid used for the discrete tapers, or check that the discrete tapers are orthogonal (the discrete tapers has a max-cross L2 norm of $max_discrete_cross)."
    end
end

"""
    check_tapers_for_data(data, tapers, bandwidth; wavenumber_res, tol, min_concentration)

This function checks if the tapers are suitable for the provided data.
In particular, it checks that the tapers work for the type of grids provided in the data.
This assumes that the base tapers have an L2 norm of 1 and are continuous.
We use the term `grid` here to refer to the sampling mechanism of a process.
So it also refer to continuously recorded data (which is the `NoGrid` type).

If you have a problem, you can use `taper_checks` to check the tapers directly.

Importantly, these checks are very expensive to compute.
In some cases failures in the checks may be due to the quality of the approximation of the Fourier transform of the tapers.
Sometimes this can be mitigated by increasing the wavenumber resolution of the tapers here.
However, when using interpolated tapers, you may also need to increase the resolution of the grid used to compute the Fourier transform of those tapers.

# Arguments
- `data`: The data to check the tapers for.
- `tapers`: The tapers to check.
- `bandwidth`: The bandwidth of the tapers.
- `wavenumber_res`: The wavenumber resolution to use for the tapers generated from the data grids.
- `tol`: The tolerance for the L2 norm of the tapers.
- `min_concentration`: The minimum concentration of the tapers before they are considered invalid and an error is thrown.
"""
function check_tapers_for_data(
        data,
        tapers,
        bandwidth;
        wavenumber_res = 1000,
        tol = 1e-1,
        min_concentration = 0.95
)
    normalisations, concentrations = taper_checks(
        data, tapers, bandwidth, wavenumber_res = wavenumber_res)
    concentrations_diag = [concentrations[i, i, k, l]
                           for i in axes(concentrations, 1),
    k in axes(concentrations, 3), l in axes(concentrations, 4)]
    concentrations_off_diag = [i == j ? zero(eltype(concentrations)) :
                               concentrations[i, j, k, l]
                               for
                               i in axes(concentrations, 1), j in axes(concentrations, 2),
    k in axes(concentrations, 3), l in axes(concentrations, 4)]
    @assert all(x -> abs(x - 1) < tol, stack(normalisations)) "All tapers must have an L2 norm of 1, but found largest difference of $(maximum(x->abs(x-1), stack(normalisations)))"
    @assert all(x -> real(x) > min_concentration, concentrations_diag) "The worst taper concentration was $(minimum(x->real(x), concentrations_diag))"
    @assert all(x -> abs(x) < tol, concentrations_off_diag) "All tapers must have no cross-concentration, but found largest difference of $(maximum(x->abs(x), concentrations_off_diag))"

    @info "Tapers are suitable for the data."
end

function taper_checks(data, tapers, bandwidth; wavenumber_res = 500)
    @assert bandwidth>0 "Bandwidth must be positive"

    observational_types = unique(get_grid.(data))
    tapers_on_grids = [tapers_on_grid(tapers, grid; wavenumber_res = wavenumber_res)
                       for grid in observational_types]
    concentration_region = Ball(Point(ntuple(i -> 0.0, _getdims(data[1]))), bandwidth)
    conc_res = choose_concentration_resolution(tapers_on_grids, concentration_region)
    normalisations = tapers_normalisations(tapers_on_grids) # checks on each family
    concentrations = taper_concentrations(
        tapers_on_grids, concentration_region; resolution = conc_res) # has to check across the sampled and non sampled versions
    return normalisations, concentrations
end

_getdims(x) = x isa PointSet ? embeddim(x) : embeddim(domain(x))
