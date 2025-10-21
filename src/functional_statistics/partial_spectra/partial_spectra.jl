"""
    partial_spectra(spectrum::Spectra{MarginalTrait}) -> Spectra{PartialTrait}
    partial_spectra(spectrum::RotationalSpectra{MarginalTrait}) -> RotationalEstimate{PartialTrait}
    partial_spectra(data, region; kwargs...) -> Spectra{PartialTrait}
    partial_spectra(data::SpatialData; kwargs...) -> Spectra{PartialTrait}

Compute partial spectral estimates from marginal spectral estimates or directly from data.

Partial spectra remove the linear influence of all other processes, providing a measure of
the direct relationship between process pairs. Unlike marginal spectra which show total
power including indirect effects, partial spectra reveal only the direct linear
relationships after removing the influence of all other processes. The computation involves
matrix inversion of the spectral matrix and finite-sample bias correction based on the
number of tapers used in the original estimation.

# Arguments
- `spectrum::Spectra{MarginalTrait}`: A marginal spectral estimate to convert
- `spectrum::RotationalSpectra{MarginalTrait}`: A rotational marginal spectral estimate
- `data`: Spatial data for direct partial spectra computation
- `region::Meshes.Geometry`: Spatial region for direct computation

# Keywords
When computing directly from data, all keywords from [`spectra`](@ref) are supported:
- `nk`: Number of wavenumbers in each dimension
- `kmax`: Maximum wavenumber in each dimension
- `dk`: Wavenumber spacing in each dimension
- `tapers`: Taper functions to use
- `nw = 3`: Space-bandwidth product for taper generation
- `mean_method = DefaultMean()`: Method for mean estimation

# Returns
- `Spectra{PartialTrait}`: Partial spectral estimate with same wavenumber grid as input
- `RotationalEstimate{PartialTrait}`: For rotational input spectra

The returned partial spectra have the same spatial dimensions and wavenumber grid as the
input, but represent direct relationships between processes rather than total power.

# Throws
- `ArgumentError`: If the spectrum does not have equal input and output process sets
    (i.e., the spectral matrix is not square). Partial spectra require square spectral
    matrices, so you shouldn't have subsetted before calling partial spectra.

# Mathematical Details
For a spectral matrix `f`, the partial spectrum `P` is computed through matrix inversion:

**Basic Formula:**
- Diagonal elements: `Pᵢᵢ = 1/Gᵢᵢ`
- Off-diagonal elements: `Pᵢⱼ = -Gᵢⱼ/(GᵢᵢGⱼⱼ - |Gᵢⱼ|²)`

where `G = f⁻¹` is the inverse spectral matrix.

**Bias Correction:**
Finite-sample bias correction is automatically applied using the number of tapers:
- Correction factor: `M/(M - Q + xᵢⱼ)` where:
  -`M = number of tapers from the original spectral estimate
  - `Q` = number of processes
  - `xᵢⱼ` = 1 if i=j (diagonal), 2 if i≠j (off-diagonal)

**Special Cases:**
- Single process: Returns the original spectral value unchanged
- Rotational spectra: Bias correction is handled differently (not fully implemented)

# Notes
- Partial spectra are the spectral domain equivalent of partial correlation
- Values can be complex and may have larger magnitudes than marginal spectra
- Diagonal elements represent the "partial power" of each process
- Off-diagonal elements show direct cross-relationships
- Use [`partial_spectra_uncorrected`](@ref) to skip bias correction
- For rotational estimates, bias correction is not currently fully implemented

# Examples
```julia
# Compute partial spectra from existing marginal estimates
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
partial_spec = partial_spectra(marginal_spec)

# Direct computation from data and region
partial_spec = partial_spectra(data, region; nk = (32, 32), kmax = (0.5, 0.5), nw = 4)

# Direct computation from SpatialData object
spatial_data_obj = spatial_data(data, region)
partial_spec = partial_spectra(spatial_data_obj; kmax = 0.3, tapers = my_tapers)

# Access partial relationships between processes
direct_coupling = partial_spec[1, 2]  # Direct coupling between processes 1 and 2
partial_power = partial_spec[1, 1]    # Partial power of process 1
```

See also: [`spectra`](@ref), [`partial_coherence`](@ref), [`partial_spectra_uncorrected`](@ref)
"""
partial_spectra

function partial_spectra(spectrum::NormalOrRotationalSpectra{MarginalTrait})
    mem = deepcopy(spectrum)
    partial_spectra!(mem)
end
function partial_spectra!(spectrum::Spectra{MarginalTrait})::Spectra{PartialTrait}
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Partial spectra computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(spectrum)). " *
            "Partial spectra measure direct relationships between process pairs."
        ))
    end

    trait = process_trait(spectrum)
    est = getestimates(spectrum)
    process_info = getprocessinformation(spectrum)
    estimation_info = getestimationinformation(spectrum)

    transformed = apply_transform!(
        _partial_spectra_noalloc!, est, trait, estimation_info.ntapers)

    return Spectra{PartialTrait}(
        getevaluationpoints(spectrum), transformed, process_info, estimation_info)
end

function partial_spectra(data, region::Meshes.Geometry; kwargs...)::Spectra{PartialTrait}
    return partial_spectra(spatial_data(data, region); kwargs...)
end

function partial_spectra(data::SpatialData; kwargs...)::Spectra{PartialTrait}
    return partial_spectra!(spectra(data; kwargs...))
end

function partial_spectra!(spectrum::RotationalSpectra{MarginalTrait})::RotationalEstimate{PartialTrait}
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Partial spectra computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(spectrum))."
        ))
    end

    wavenumber = getevaluationpoints(spectrum)
    power = getestimates(spectrum)
    trait = process_trait(spectrum)
    # Note: Debiasing is different for rotational estimates (not currently implemented)
    value = apply_transform!(_partial_spectra_noalloc!, power, trait, nothing)
    processinfo = getprocessinformation(spectrum)
    estimationinfo = getestimationinformation(spectrum)

    return RotationalEstimate{PartialTrait, typeof(spectrum)}(
        wavenumber, value, processinfo, estimationinfo)
end

function partial_spectra(x::SMatrix, ::Nothing)
    g = inv(x)
    A = diagm(diag(g))
    g2 = abs2.(g)
    denom = A * ones(typeof(x)) * A - g2 + diagm(diag(g2))
    # Compute partial spectra: -gⱼₖ / (gⱼⱼ gₖₖ - |gⱼₖ|²) if j ≠ k, 1 / gⱼⱼ if j = k
    return (g ./ denom) .* (2I - ones(typeof(x)))
end

function partial_spectra!(x::AbstractMatrix, ::Nothing)
    C = LinearAlgebra.inv!(cholesky!(x))
    for i in axes(C, 1), j in axes(C, 2)
        if i == j
            continue
        end
        C[i, j] = -C[i, j] / (C[i, i] * C[j, j] - abs2(C[i, j]))
    end
    # need to not overwrite C[i,i] in the above loop
    for i in axes(C, 1)
        C[i, i] = 1 / C[i, i]
    end
    return C
end

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

function partial_spectra(x::SMatrix{Q, Q, T, N}, ntapers::Int) where {Q, T, N}
    p = partial_spectra(x, nothing)
    # Bias correction: M-Q+1 for diagonal, M-Q+2 for off-diagonal
    denom = ntapers .* ones(typeof(x)) .- Q .+ 2 - I
    return ntapers ./ denom .* p
end

function partial_spectra!(x::AbstractMatrix{T}, ntapers::Int) where {T}
    Q = size(x, 1)
    p = partial_spectra!(x, nothing)

    # Apply bias correction element-wise
    @inbounds for i in axes(p, 1), j in axes(p, 2)
        correction_factor = ntapers / (ntapers - Q + 2 - (i == j))
        p[i, j] *= correction_factor
    end
    return p
end

function partial_spectra(x::Number, ntapers)
    return x
end

_partial_spectra_noalloc!(x::Union{Number, SMatrix}, ntapers) = partial_spectra(x, ntapers)
_partial_spectra_noalloc!(x::AbstractMatrix, ntapers) = partial_spectra!(x, ntapers)

"""
    partial_spectra_uncorrected(spectrum::Spectra{MarginalTrait}) -> Spectra{PartialTrait}

Compute partial spectra without finite-sample bias correction.

This function computes partial spectra using the raw inverse relationship without the
finite-sample bias correction that accounts for the number of tapers. Results may be
biased for small numbers of tapers, so use with caution. Prefer [`partial_spectra`](@ref)
for most applications.

# Arguments
- `spectrum::Spectra{MarginalTrait}`: A marginal spectral estimate

# Returns
- `Spectra{PartialTrait}`: Uncorrected partial spectral estimates

# Notes
- Equivalent to calling `partial_spectra(x, nothing)` on the spectral matrices
- Useful for theoretical analysis or when bias correction is undesired
- Results will differ from corrected partial spectra, especially with few tapers

# Examples
```julia
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
uncorrected_partial = partial_spectra_uncorrected(marginal_spec)
corrected_partial = partial_spectra(marginal_spec)
```

See also: [`partial_spectra`](@ref)
"""
function partial_spectra_uncorrected(spectrum::Spectra{MarginalTrait})::Spectra{PartialTrait}
    # Create a modified spectrum with no taper information for uncorrected computation
    new_spectrum = Spectra{MarginalTrait}(
        getevaluationpoints(spectrum), getestimates(spectrum),
        getprocessinformation(spectrum), EstimationInformation(nothing))
    return partial_spectra(new_spectrum)
end
