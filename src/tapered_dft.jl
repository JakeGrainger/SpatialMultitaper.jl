"""
    tapered_dft(data::SpatialData, tapers, nk, kmax, mean_method=DefaultMean())

Compute the tapered discrete Fourier transform of spatial data.

Returns different output formats depending on the input data structure:

# Single SpatialData
For a single process, returns a single array of the transformed data.

# MultipleSpatialDataVec
For `P` processes stored as a vector, returns an array of size
`P × M × n1 × ... × nD`, where `M` is the number of tapers and
`n1, ..., nD` are the wavenumber dimensions.

# MultipleSpatialDataTuple
For `P` processes stored as a tuple, returns a tuple of length `P`, where each
entry is an array of size `n1 × ... × nD × M`.

# Arguments
- `data::SpatialData`: The spatial data to transform
- `tapers`: Collection of taper functions to apply
- `nk`: Number of wavenumbers in each dimension
- `kmax`: Maximum wavenumber in each dimension
- `mean_method::MeanEstimationMethod = DefaultMean()`: Method for mean estimation

# Notes
The `mean_method` should have the same length as the number of processes in `data`.
"""
function tapered_dft(
        data, tapers, nk, kmax, mean_method::MeanEstimationMethod = DefaultMean())
    mem = preallocate_tapered_dft(data, tapers, nk, kmax)
    return tapered_dft!(mem, data, tapers, nk, kmax, mean_method)
end

function tapered_dft!(mem, sd::SingleProcessData, tapers, nk, kmax,
        mean_method::MeanEstimationMethod = DefaultMean())
    checkmeanmethod(sd, mean_method)
    return _single_tapered_dft!(mem[1], mem[2], sd, tapers, nk, kmax, mean_method)
end

function tapered_dft!(mem, sd::MultipleSpatialDataVec, tapers, nk, kmax,
        mean_method::MeanEstimationMethod = DefaultMean())
    checkmeanmethod(sd, mean_method)
    for p in 1:ncol(sd)
        result = _single_tapered_dft!(
            mem[1][p], mem[2][p], sd[p], tapers, nk, kmax, mean_method[p])
        for idx in CartesianIndices(size(result)[1:(end - 1)])
            for m in axes(result, ndims(result))
                mem[3][p, m, idx] = result[idx, m]
            end
        end
    end
    return mem[3]
end

function tapered_dft!(mem, sd::MultipleSpatialDataTuple{P}, tapers, nk, kmax,
        mean_method::MeanEstimationMethod = DefaultMean()) where {P}
    checkmeanmethod(sd, mean_method)
    dfts = ntuple(
        p -> _single_tapered_dft!(
            mem[1][p], mem[2][p], sd[p], tapers, nk, kmax, mean_method[p]),
        Val{P}())
    return dfts
end

function preallocate_tapered_dft(data::SingleProcessData, tapers, nk, kmax)
    tapered_data = preallocate_single_tapered_dft(data, tapers)
    mem = preallocate_fft_any_type(tapered_data, data, nk, kmax)
    return (mem, tapered_data)
end

function preallocate_tapered_dft(data::MultipleSpatialDataVec, tapers, nk, kmax)
    tapered_data = [preallocate_single_tapered_dft(data[p], tapers) for p in 1:ncol(data)]
    mem = [preallocate_fft_any_type(tapered_data[p], data[p], nk, kmax)
           for p in 1:ncol(data)]
    stacked_output = zeros(complex(Float64), ncol(data), length(tapers), nk...) # TODO: should be made generic
    return (mem, tapered_data, stacked_output)
end

function preallocate_tapered_dft(
        data::MultipleSpatialDataTuple{P}, tapers, nk, kmax) where {P}
    tapered_data = ntuple(
        p -> preallocate_single_tapered_dft(data[p], tapers), Val{P}())
    mem = ntuple(
        p -> preallocate_fft_any_type(tapered_data[p], data[p], nk, kmax), Val{P}())
    return (mem, tapered_data)
end

##

function _single_tapered_dft(data, tapers, nk, kmax, mean_method::MeanEstimationMethod)
    tapered_data = preallocate_single_tapered_dft(data, tapers)
    mem = preallocate_fft_any_type(tapered_data, data, nk, kmax)
    return _single_tapered_dft!(mem, tapered_data, data, tapers, nk, kmax, mean_method)
end

function preallocate_fft_any_type(tapered_data, data::PointPattern, nk, kmax)
    points = observations(data)
    mem = precompute_nufft_anydomain_output(getregion(data), nk, kmax, points, tapered_data)
    return mem
end

function preallocate_fft_any_type(tapered_data, data::MarkedPointPattern, nk, kmax)
    points = domain(observations(data))
    mem = precompute_nufft_anydomain_output(getregion(data), nk, kmax, points, tapered_data)
    return mem
end

function preallocate_fft_any_type(tapered_data, data::GriddedData, nk, kmax)
    grid = domain(observations(data))
    mem = preallocate_fft_anydomain(tapered_data, grid, nk, kmax)
    return mem
end

function preallocate_single_tapered_dft(data::PointPattern, tapers)
    n = length(observations(data))
    marks = one(Float64) # TODO: should be made generic
    return Vector{complex(eltype(marks))}(undef, n * length(tapers))
end

function preallocate_single_tapered_dft(data::MarkedPointPattern, tapers)
    n = length(domain(observations(data)))
    marks = values(observations(data))[1]
    return Vector{complex(eltype(marks))}(undef, n * length(tapers))
end

function _single_tapered_dft!(
        mem, tapered_marks, data::Union{PointPattern, MarkedPointPattern}, tapers,
        nk, kmax, mean_method::MeanEstimationMethod)
    # apply tapers
    λ = mean_estimate(data, mean_method)

    points, marks = _unpack_points_and_marks(data)
    region = getregion(data)
    wavenumber = Iterators.ProductIterator(_make_wavenumber_grid(nk, kmax))

    _apply_taper!(tapered_marks, points, marks, tapers)

    # perform transform
    J = nufft_anydomain!(mem, region, nk, kmax, points, tapered_marks, -1, 1e-14)
    for i in eachindex(tapers)
        J_h = selectdim(J, embeddim(points) + 1, i)
        for (j, k) in zip(eachindex(J_h), wavenumber)
            J_h[j] -= λ * taper_ft(tapers[i], k)
        end
    end
    return J
end
function _unpack_points_and_marks(data::MarkedPointPattern)
    marks = values(observations(data))[1]
    points = domain(observations(data))
    return points, marks
end
function _unpack_points_and_marks(data::PointPattern)
    marks = one(Float64) # TODO: should be made generic
    points = observations(data)
    return points, marks
end

function _apply_taper!(tapered_marks, points::PointSet, marks, tapers)
    n = length(points)
    for i in eachindex(tapers)
        tapered_marks[(1:n) .+ (i - 1) * n] .= complex.(tapers[i].(points) .* marks)
    end
    return tapered_marks
end

function preallocate_single_tapered_dft(data::GriddedData, tapers)
    rf = values(observations(data))[1]
    grid = domain(observations(data))
    return zeros(eltype(rf), (size(grid)..., length(tapers)))
end

function _single_tapered_dft!(mem, tapered_data, data::GriddedData, tapers,
        nk, kmax, mean_method::MeanEstimationMethod)
    λ = mean_estimate(data, mean_method)
    rf = values(observations(data))[1]
    grid = domain(observations(data))

    # apply taper and mean removal
    scaling = prod(unitless_spacing(grid))
    for i in eachindex(tapers)
        for (j, idx) in enumerate(CartesianIndices(size(grid)))
            rf_val = rf[j]
            rf_val = isnan(rf_val) ? zero(eltype(rf)) : rf_val
            scaled_rf = (rf_val - λ) * scaling
            loc = centroid(grid, j)
            val = tapers[i](loc) * scaled_rf
            tapered_data[idx, i] = val
        end
    end

    # perform transform
    return fft_anydomain!(mem, tapered_data, grid, nk, kmax)
end
