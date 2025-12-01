computed_from(::Type{<:Spectra{PartialTrait, D}}) where {D} = Spectra{MarginalTrait, D}
function computed_from(::Type{<:RotationalSpectra{PartialTrait, D}}) where {D}
    Spectra{PartialTrait, D}
end

function allocate_estimate_memory(::Type{<:NormalOrRotationalSpectra{PartialTrait}},
        ::Type{<:NormalOrRotationalSpectra{MarginalTrait}}, relevant_memory; kwargs...)
    return relevant_memory, nothing
end

function extract_relevant_memory(::Type{<:NormalOrRotationalSpectra{PartialTrait}},
        est::NormalOrRotationalSpectra{MarginalTrait})
    return deepcopy(get_estimates(est))
end
function extract_relevant_memory(::Type{<:NormalOrRotationalSpectra{PartialTrait}},
        est::EstimateMemory{<:NormalOrRotationalSpectra{MarginalTrait}})
    return est.output_memory
end

function validate_core_parameters(
        ::Type{<:Spectra{PartialTrait}}; kwargs...)
    # no additional parameters to validate
    return nothing
end

function resolve_missing_parameters(
        ::Type{<:Spectra{PartialTrait}}, arg; kwargs...)
    return kwargs
end

function validate_memory_compatibility(
        ::Type{<:Spectra{PartialTrait}}, mem, source; kwargs...)
    # no additional compatibility checks needed
    @argcheck size(mem.output_memory) == size(get_estimates(source))
    @argcheck eltype(mem.output_memory) == eltype(get_estimates(source))
end

function compute_estimate!(::Type{<:NormalOrRotationalSpectra{PartialTrait}},
        mem, source::NormalOrRotationalSpectra{MarginalTrait};
        kwargs...)
    if !is_same_process_sets(source)
        throw(ArgumentError(
            "Partial spectra computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(source))."
        ))
    end

    trait = process_trait(source)
    est = get_estimates(source)
    process_info = get_process_information(source)
    estimation_info = get_estimation_information(source)

    transformed = apply_transform!(
        _partial_spectra_noalloc!, mem.output_memory, est, trait, estimation_info.ntapers)

    return _build_partial_output(
        source, get_evaluation_points(source), transformed, process_info, estimation_info)
end

function _build_partial_output(::Spectra{MarginalTrait}, args...)
    return Spectra{PartialTrait}(args...)
end
function _build_partial_output(
        ::RotationalSpectra{MarginalTrait, D, S}, args...) where {D, S}
    return RotationalSpectra{PartialTrait, S}(args...)
end

## internals

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
