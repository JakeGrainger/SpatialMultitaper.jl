##
using SpatialMultitaper

function _internal_c(out, store, data; radii, kwargs...)
    spectrum = partial_spectra(data; kwargs...)
    value = SpatialMultitaper._sdf2C!(out, store, spectrum, radii)
end

function foo()
    points1 = PointSet(Point.(rand(1000), rand(1000)))
    points2 = PointSet(Point.(rand(1000), rand(1000)))
    points3 = PointSet(Point.(rand(1000), rand(1000)))
    box = Box((0, 0), (1, 1))
    data = spatial_data((points1, points2, points3), box)

    grid1 = georef(
        (rf = rand(100 * 100),), CartesianGrid((0, 0), (1, 1), dims = (100, 100)))
    grid2 = georef(
        (rf = rand(100 * 100),), CartesianGrid((0, 0), (1, 1), dims = (100, 100)))
    grid3 = georef(
        (rf = rand(100 * 100),), CartesianGrid((0, 0), (1, 1), dims = (100, 100)))
    data = spatial_data((grid1, grid2, grid3), box)

    radii = range(0, 0.25, 100)
    tapers = sin_taper_family((4, 4), box)
    spectrum = partial_spectra(data, kmax = 100, nk = 100, tapers = tapers)
    out = SpatialMultitaper.preallocate_c_output(spectrum, radii)
    store = SpatialMultitaper.precompute_c_weights(spectrum, radii)
    @profview _internal_c(
        out, store, data; radii, kmax = 100, nk = 100, tapers = tapers)
    #  partial_l_function(data, radii = radii, kmax = 100, nk = 100, tapers = tapers)
end

foo()
