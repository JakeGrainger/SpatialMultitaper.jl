"""
	tapered_dft(data::AbstractVector{<:Union{GeoTable,PointSet}}, tapers, nfreq, fmax, region, mean_method::AbstractVector{<:MeanEstimationMethod})

Compute the tapered discrete Fourier transform of a collection of data sets.
If we have `n_1,...,n_D` frequencies, and M tapers and P processes, this returns an array of size `P x M x n_1 x ... x n_D`.
Note that for small numbers of processes, it is more efficient to use the `NTuple` method.
"""
function tapered_dft(
    data::AbstractVector{<:Union{GeoTable,PointSet}},
    tapers,
    nfreq,
    fmax,
    region,
    mean_method::AbstractVector{<:MeanEstimationMethod},
)
    @assert length(data) == length(mean_method) "data and mean_method must have the same length, got $(length(data)) and $(length(mean_method))."
    dfts = stack(
        single_tapered_dft(data[p], tapers, nfreq, fmax, region, mean_method[p]) for
        p in eachindex(data)
    )
    return permutedims(dfts, (ndims(dfts), ndims(dfts) - 1, 1:ndims(dfts)-2...))
end

"""
	tapered_dft(data::NTuple{P,Union{GeoTable, PointSet}}, tapers, nfreq, fmax, region, mean_method::NTuple{P,MeanEstimationMethod}) where {P}

Compute the tapered discrete Fourier transform of a collection of data sets.
If we have n_1,...,n_D frequencies, and M tapers and P processes, this returns a `Tuple` of size `P`, whose entries are arrays of size `n_1 x ... x n_D x M`.
"""
function tapered_dft(
    data::NTuple{P,Union{GeoTable,PointSet}},
    tapers,
    nfreq,
    fmax,
    region,
    mean_method::NTuple{P,MeanEstimationMethod},
) where {P}
    dfts = ntuple(
        p -> single_tapered_dft(data[p], tapers, nfreq, fmax, region, mean_method[p]),
        Val{P}(),
    )
    return dfts
end

# points (just adds marks with value 1 and proceeds)
function single_tapered_dft(
    data::PointSet,
    tapers,
    nfreq,
    fmax,
    region,
    mean_method::MeanEstimationMethod,
)
    marks = [1.0 * (data[i] ∈ region) for i ∈ eachindex(data)]
    single_tapered_dft(data, marks, tapers, nfreq, fmax, region, mean_method)
end

# allows for dispatch with new GeoTable
function single_tapered_dft(
    data,
    tapers,
    nfreq,
    fmax,
    region,
    mean_method::MeanEstimationMethod,
)
    single_tapered_dft(
        domain(data),
        values(data)[1],
        tapers,
        nfreq,
        fmax,
        region,
        mean_method,
    )
end

# points with marks
function single_tapered_dft(
    points::PointSet,
    marks,
    tapers,
    nfreq,
    fmax,
    region,
    mean_method::MeanEstimationMethod,
)
    # apply tapers
    freq = Iterators.ProductIterator(make_freq(nfreq, fmax, embeddim(points)))
    tapered_marks = apply_taper(points, marks, tapers)

    # perform transform
    λ = mean_estimate(points, marks, region, mean_method)
    J = nufft_anydomain(region, nfreq, fmax, points, tapered_marks, -1, 1e-14)
    for i ∈ eachindex(tapers)
        J_h = selectdim(J, embeddim(points) + 1, i)
        for (j, k) ∈ zip(eachindex(J_h), freq)
            J_h[j] -= λ * taper_ft(tapers[i], k)
        end
    end
    return J
end
function apply_taper(points::PointSet, marks, tapers)
    n = length(points)
    tapered_marks = Vector{complex(eltype(marks))}(undef, n * length(tapers))
    for i ∈ eachindex(tapers)
        tapered_marks[(1:n).+(i-1)*n] .= complex.(tapers[i].(points) .* marks)
    end
    return tapered_marks
end

function single_tapered_dft(
    grid::CartesianGrid,
    rf,
    tapers,
    nfreq,
    fmax,
    region,
    mean_method::MeanEstimationMethod,
)
    # preallocate    
    tapered_data =
        Array{eltype(rf),embeddim(grid) + 1}(undef, (size(grid)..., length(tapers)))
    grid_points = Iterators.ProductIterator(grid2side(grid))

    # apply taper and mean removal
    λ = mean_estimate(grid, rf, region, mean_method)
    scaling = prod(unitless_spacing(grid))
    for i ∈ eachindex(tapers)
        for (j, u) ∈
            zip(eachindex(selectdim(tapered_data, embeddim(grid) + 1, i)), grid_points)
            selectdim(tapered_data, embeddim(grid) + 1, i)[j] =
                tapers[i](u) * (rf[j] - λ) * scaling
        end
    end

    # perform transform
    J = fft_anydomain(tapered_data, grid, nfreq, fmax)
    return J
end

