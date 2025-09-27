function partial_spectra(spectrum::Spectra{MarginalTrait})
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Partial spectra can only be computed for estimates where the two sets of " *
            "processes are the same."
        ))
    end
    trait = estimate_trait(spectrum)
    est = getestimate(spectrum)
    process_info = getprocessinformation(spectrum)
    estimation_info = getestimationinformation(spectrum)
    transformed = apply_transform(partial_spectra, est, trait, estimation_info.ntapers)
    return Spectra{PartialTrait}(
        getargument(spectrum), transformed, process_info, estimation_info)
end

function partial_spectra(data, region; kwargs...)
    partial_spectra(spatial_data(data, region); kwargs...)
end
partial_spectra(data::SpatialData; kwargs...) = partial_spectra(spectra(data; kwargs...))

function partial_spectra_uncorrected(spectrum::Spectra{MarginalTrait})
    new_spectrum = Spectra{MarginalTrait}(getargument(spectrum), getestimate(spectrum),
        getprocessinformation(spectrum), EstimationInformation(nothing))
    partial_spectra(new_spectrum)
end

function partial_spectra(spectrum::RotationalSpectra{MarginalTrait})
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Partial spectra can only be computed for estimates where the two sets of " *
            "processes are the same."
        ))
    end
    freq = getargument(spectrum)
    power = getestimate(spectrum)
    trait = estimate_trait(spectrum)
    value = apply_transform(partial_spectra, power, trait, nothing) # debiasing is different here (not currently implemented)
    processinfo = getprocessinformation(spectrum)
    estimationinfo = getestimationinformation(spectrum)
    return RotationalEstimate{PartialTrait, typeof(spectrum)}(
        freq, value, processinfo, estimationinfo)
end

function partial_spectra(x::SMatrix, ::Nothing)
    g = inv(x)
    A = diagm((diag(g)))
    g2 = abs2.(g)
    denom = A * ones(typeof(x)) * A - g2 + diagm(diag(g2))
    return (g ./ denom) .* (2I - ones(typeof(x)))
    # computes -gⱼₖ / (gⱼⱼ gₖₖ - |gⱼₖ|²) if j ≠ k
    # computes  1 / gⱼⱼ if j = k
end

function partial_spectra(x::AbstractMatrix, ::Nothing)
    C = inv(x)

    return [i == j ? 1 / C[i, i] : -C[i, j] / (C[i, i] * C[j, j] - abs2(C[i, j]))
            for
            i in axes(C, 1), j in axes(C, 2)]
end

function partial_spectra(x::SMatrix{2, 2, T, 4}, ::Nothing) where {T}
    g = inv(x)
    p11 = 1 / g[1, 1]
    p22 = 1 / g[2, 2]
    p12 = -g[1, 2] / (g[1, 1] * g[2, 2] - abs2(g[1, 2]))
    p21 = conj(p12)

    return SMatrix{2, 2, T, 4}(
        p11,
        p21,
        p12,
        p22 # column major
    )
end

function partial_spectra(x::SMatrix{Q, Q, T, N}, ntapers::Int) where {Q, T, N}
    p = partial_spectra(x, nothing)
    denom = ntapers .* ones(typeof(x)) .- Q .+ 2 - I # so that M - Q + 2 off diag and M - Q + 1 on diag
    return ntapers ./ denom .* p
end

function partial_spectra(x::AbstractMatrix{T}, ntapers::Int) where {T}
    Q = size(x, 1)
    p = partial_spectra(x, nothing)
    for i in axes(p, 1), j in axes(p, 2)
        @inbounds p[i, j] *= ntapers / (ntapers - Q + 2 - (i == j))
    end
    return p
end

function partial_spectra(x::Number, ntapers)
    return x
end
