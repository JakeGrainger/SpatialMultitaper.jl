"""
    tapered_dft(data::SpatialData, tapers, nfreq, fmax, mean_method=DefaultMean())

Compute the tapered discrete Fourier transform of spatial data.

Returns different output formats depending on the input data structure:

# Single SpatialData
For a single process, returns a single array of the transformed data.

# MultipleSpatialDataVec
For `P` processes stored as a vector, returns an array of size
`P × M × n1 × ... × nD`, where `M` is the number of tapers and
`n1, ..., nD` are the frequency dimensions.

# MultipleSpatialDataTuple
For `P` processes stored as a tuple, returns a tuple of length `P`, where each
entry is an array of size `n1 × ... × nD × M`.

# Arguments
- `data::SpatialData`: The spatial data to transform
- `tapers`: Collection of taper functions to apply
- `nfreq`: Number of frequencies in each dimension
- `fmax`: Maximum frequency in each dimension
- `mean_method::MeanEstimationMethod = DefaultMean()`: Method for mean estimation

# Notes
The `mean_method` should have the same length as the number of processes in `data`.
"""
function tapered_dft(sd::SpatialData, tapers, nfreq, fmax,
        mean_method::MeanEstimationMethod = DefaultMean())
    checkmeanmethod(sd, mean_method)
    return _single_tapered_dft(sd, tapers, nfreq, fmax, mean_method)
end
function tapered_dft(sd::MultipleSpatialDataVec, tapers, nfreq, fmax,
        mean_method::MeanEstimationMethod = DefaultMean())
    checkmeanmethod(sd, mean_method)
    dfts = stack(_single_tapered_dft(sd[p], tapers, nfreq, fmax, mean_method[p])
    for p in 1:ncol(sd))
    return permutedims(dfts, (ndims(dfts), ndims(dfts) - 1, 1:(ndims(dfts) - 2)...))
end
function tapered_dft(sd::MultipleSpatialDataTuple, tapers, nfreq, fmax,
        mean_method::MeanEstimationMethod = DefaultMean())
    checkmeanmethod(sd, mean_method)
    dfts = ntuple(
        p -> _single_tapered_dft(sd[p], tapers, nfreq, fmax, mean_method[p]),
        Val{ncol(sd)}()
    )
    return dfts
end

# points (just adds marks with value 1 and proceeds)
function _single_tapered_dft(
        data::PointPattern, tapers, nfreq, fmax, mean_method::MeanEstimationMethod)
    marks = ones(length(observations(data)))
    markedtable = georef((marks = marks,), observations(data))
    markeddata = SpatialData(markedtable, getregion(data), propertynames(data))
    return _single_tapered_dft(markeddata, tapers, nfreq, fmax, mean_method)
end

# points with marks
function _single_tapered_dft(
        data::MarkedPointPattern, tapers, nfreq, fmax, mean_method::MeanEstimationMethod)
    # apply tapers
    λ = mean_estimate(data, mean_method)

    marks = values(observations(data))[1]
    points = domain(observations(data))
    region = getregion(data)
    freq = Iterators.ProductIterator(make_freq(nfreq, fmax, embeddim(points)))
    tapered_marks = _apply_taper(points, marks, tapers)

    # perform transform
    J = nufft_anydomain(region, nfreq, fmax, points, tapered_marks, -1, 1e-14)
    for i in eachindex(tapers)
        J_h = selectdim(J, embeddim(points) + 1, i)
        for (j, k) in zip(eachindex(J_h), freq)
            J_h[j] -= λ * taper_ft(tapers[i], k)
        end
    end
    return J
end
function _apply_taper(points::PointSet, marks, tapers)
    n = length(points)
    tapered_marks = Vector{complex(eltype(marks))}(undef, n * length(tapers))
    for i in eachindex(tapers)
        tapered_marks[(1:n) .+ (i - 1) * n] .= complex.(tapers[i].(points) .* marks)
    end
    return tapered_marks
end

function _single_tapered_dft(data::GriddedData, tapers, nfreq,
        fmax, mean_method::MeanEstimationMethod)
    # preallocate
    λ = mean_estimate(data, mean_method)
    rf = values(observations(data))[1]
    grid = domain(observations(data))

    replace!(rf, NaN => zero(eltype(rf))) # means single_tapered_dft should have ! after name

    tapered_data = Array{eltype(rf), embeddim(grid) + 1}(
        undef, (size(grid)..., length(tapers)))

    # apply taper and mean removal
    scaling = prod(unitless_spacing(grid))
    for i in eachindex(tapers)
        for j in eachindex(selectdim(tapered_data, embeddim(grid) + 1, i))
            val = tapers[i](centroid(grid, j)) * (rf[j] - λ) * scaling
            selectdim(tapered_data, embeddim(grid) + 1, i)[j] = val
        end
    end

    # perform transform
    J = fft_anydomain(tapered_data, grid, nfreq, fmax)
    return J
end
