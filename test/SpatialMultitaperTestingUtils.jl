module SpatialMultitaperTestingUtils
using SpatialMultitaper
function slow_dft(u, f, freq, iflag)
    pm = iflag ≥ 0 ? 1 : -1
    return [
        sum(f[i] * exp(pm * 2pi * 1im * sum(u[i] .* k)) for i in eachindex(u, f)) for
        k in freq
    ]
end

function make_simple_example(; vector_of_processes = false)
    region = Box(Point(0, 0), Point(3, 3))
    pattern = georef((marks = [1, 2.4],), PointSet([Point(0, 0), Point(1, 1)]))
    pattern2 = PointSet([Point(0, 0.3), Point(0.2, 0.2), Point(0.1, 1)])
    griddata = georef((rf = rand(6 * 6),), CartesianGrid((0, 0), (3, 3), dims = (6, 6)))
    tapers = sin_taper_family((3, 3), region)
    nfreq = (10, 10)
    fmax = (2, 2)
    data =
        vector_of_processes ? [pattern, pattern2, griddata] : (pattern, pattern2, griddata)
    mt_est = multitaper_estimate(data, region, nfreq = nfreq, fmax = fmax, tapers = tapers)
    return mt_est
end
export slow_dft, make_simple_example
end