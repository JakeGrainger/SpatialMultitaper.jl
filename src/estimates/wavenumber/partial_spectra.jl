"""
    partial_spectra(spectrum::Spectra{MarginalTrait})

Compute partial spectral estimates from marginal spectral estimates.

Partial spectra remove the linear influence of all other processes, providing
a measure of the direct relationship between process pairs. The computation
involves matrix inversion and bias correction based on the number of tapers.

# Arguments
- `spectrum::Spectra{MarginalTrait}`: A marginal spectral estimate

# Returns
A `Spectra{PartialTrait}` object containing the partial spectral estimates.

# Throws
- `ArgumentError`: If the spectrum does not have equal input and output process sets

# Mathematical Details
For a spectral matrix S, the partial spectrum P is computed as:
- Pᵢᵢ = 1/Gᵢᵢ (diagonal elements)
- Pᵢⱼ = -Gᵢⱼ/(GᵢᵢGⱼⱼ - |Gᵢⱼ|²) (off-diagonal elements)
where G = S⁻¹ is the inverse spectral matrix.

# Examples
```julia
# Compute partial spectra from marginal estimates
partial_spec = partial_spectra(marginal_spec)

# Direct computation from data::SpatialData
partial_spec = partial_spectra(data; nk=nk, kmax=kmax, tapers=tapers)
```
"""
function partial_spectra(spectrum::Spectra{MarginalTrait})::Spectra{PartialTrait}
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Partial spectra computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(spectrum)). " *
            "Partial spectra measure direct relationships between process pairs."
        ))
    end

    trait = process_trait(spectrum)
    est = getestimate(spectrum)
    process_info = getprocessinformation(spectrum)
    estimation_info = getestimationinformation(spectrum)

    transformed = apply_transform(partial_spectra, est, trait, estimation_info.ntapers)

    return Spectra{PartialTrait}(
        getargument(spectrum), transformed, process_info, estimation_info)
end

"""
    partial_spectra(data, region; kwargs...)

Compute partial spectral estimates directly from data and region.
"""
function partial_spectra(data, region; kwargs...)::Spectra{PartialTrait}
    return partial_spectra(spatial_data(data, region); kwargs...)
end

"""
    partial_spectra(data::SpatialData; kwargs...)

Compute partial spectral estimates from spatial data.

First computes marginal spectra, then converts to partial spectra.
"""
function partial_spectra(data::SpatialData; kwargs...)::Spectra{PartialTrait}
    return partial_spectra(spectra(data; kwargs...))
end

"""
    partial_spectra_uncorrected(spectrum::Spectra{MarginalTrait})

Compute partial spectra without bias correction from the number of tapers.

This function computes partial spectra using the raw inverse relationship
without the finite-sample bias correction that accounts for the number of tapers.
Use with caution as results may be biased for small numbers of tapers.

# Arguments
- `spectrum::Spectra{MarginalTrait}`: A marginal spectral estimate

# Returns
A `Spectra{PartialTrait}` object with uncorrected partial spectral estimates.
"""
function partial_spectra_uncorrected(spectrum::Spectra{MarginalTrait})::Spectra{PartialTrait}
    # Create a modified spectrum with no taper information for uncorrected computation
    new_spectrum = Spectra{MarginalTrait}(getargument(spectrum), getestimate(spectrum),
        getprocessinformation(spectrum), EstimationInformation(nothing))
    return partial_spectra(new_spectrum)
end

"""
    partial_spectra(spectrum::RotationalSpectra{MarginalTrait})

Compute partial spectral estimates from rotational marginal spectral estimates.

For rotational estimates, the bias correction is handled differently and is
not currently fully implemented.
"""
function partial_spectra(spectrum::RotationalSpectra{MarginalTrait})::RotationalEstimate{PartialTrait}
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Partial spectra computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(spectrum))."
        ))
    end

    freq = getargument(spectrum)
    power = getestimate(spectrum)
    trait = process_trait(spectrum)
    # Note: Debiasing is different for rotational estimates (not currently implemented)
    value = apply_transform(partial_spectra, power, trait, nothing)
    processinfo = getprocessinformation(spectrum)
    estimationinfo = getestimationinformation(spectrum)

    return RotationalEstimate{PartialTrait, typeof(spectrum)}(
        freq, value, processinfo, estimationinfo)
end

"""
    partial_spectra(x::SMatrix, ::Nothing)

Compute uncorrected partial spectra from a static matrix.

This is the core mathematical transformation that converts a spectral matrix
to partial spectral matrix through matrix inversion and element-wise operations.

# Mathematical Formula
For spectral matrix S with inverse G = S⁻¹:
- Diagonal: Pᵢᵢ = 1/Gᵢᵢ
- Off-diagonal: Pᵢⱼ = -Gᵢⱼ/(GᵢᵢGⱼⱼ - |Gᵢⱼ|²)
"""
function partial_spectra(x::SMatrix, ::Nothing)
    g = inv(x)
    A = diagm(diag(g))
    g2 = abs2.(g)
    denom = A * ones(typeof(x)) * A - g2 + diagm(diag(g2))
    # Compute partial spectra: -gⱼₖ / (gⱼⱼ gₖₖ - |gⱼₖ|²) if j ≠ k, 1 / gⱼⱼ if j = k
    return (g ./ denom) .* (2I - ones(typeof(x)))
end

"""
    partial_spectra(x::AbstractMatrix, ::Nothing)

Compute uncorrected partial spectra from a general matrix.

Uses explicit indexing for maximum clarity and type stability.
"""
function partial_spectra(x::AbstractMatrix, ::Nothing)
    C = inv(x)
    out = zeros(eltype(C), size(C))
    @inbounds for i in axes(C, 1), j in axes(C, 2)
        out[i, j] = i == j ? 1 / C[i, i] : -C[i, j] / (C[i, i] * C[j, j] - abs2(C[i, j]))
    end
    return out
end

"""
    partial_spectra(x::SMatrix{2, 2, T, 4}, ::Nothing) where {T}

Optimized partial spectra computation for 2×2 static matrices.

This specialized method provides optimized computation for the common case
of two-process analysis with explicit element-wise calculations.
"""
function partial_spectra(x::SMatrix{2, 2, T, 4}, ::Nothing) where {T}
    g = inv(x)

    # Compute partial spectra elements explicitly for 2x2 case
    p11 = 1 / g[1, 1]
    p22 = 1 / g[2, 2]
    p12 = -g[1, 2] / (g[1, 1] * g[2, 2] - abs2(g[1, 2]))
    p21 = conj(p12)  # Hermitian symmetry

    return SMatrix{2, 2, T, 4}(
        p11, p21, p12, p22  # column major storage
    )
end

"""
    partial_spectra(x::SMatrix{Q, Q, T, N}, ntapers::Int) where {Q, T, N}

Compute bias-corrected partial spectra from a static matrix.

Applies finite-sample bias correction based on the number of tapers used
in the original spectral estimation.

# Arguments
- `x`: Square spectral matrix
- `ntapers`: Number of tapers used in original estimation

# Mathematical Details
The bias correction factor is (M)/(M - Q + δᵢⱼ) where:
- M = number of tapers
- Q = number of processes
- δᵢⱼ = 1 if i=j (diagonal), 2 if i≠j (off-diagonal)
"""
function partial_spectra(x::SMatrix{Q, Q, T, N}, ntapers::Int) where {Q, T, N}
    p = partial_spectra(x, nothing)
    # Bias correction: M-Q+1 for diagonal, M-Q+2 for off-diagonal
    denom = ntapers .* ones(typeof(x)) .- Q .+ 2 - I
    return ntapers ./ denom .* p
end

"""
    partial_spectra(x::AbstractMatrix{T}, ntapers::Int) where {T}

Compute bias-corrected partial spectra from a general matrix.

# Arguments
- `x`: Square spectral matrix
- `ntapers`: Number of tapers used in original estimation
"""
function partial_spectra(x::AbstractMatrix{T}, ntapers::Int) where {T}
    Q = size(x, 1)
    p = partial_spectra(x, nothing)

    # Apply bias correction element-wise
    @inbounds for i in axes(p, 1), j in axes(p, 2)
        correction_factor = ntapers / (ntapers - Q + 2 - (i == j))
        p[i, j] *= correction_factor
    end
    return p
end

"""
    partial_spectra(x::Number, ntapers)

Partial spectra for a single process (identity operation).

For a single process, the partial spectrum is identical to the original spectrum
since there are no other processes to partial out.
"""
function partial_spectra(x::Number, ntapers)
    return x
end
