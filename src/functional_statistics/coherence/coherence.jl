
"""
    Coherence{E, D, N, A, T, IP, IE} <: AbstractEstimate{E, D, N}

A coherence estimate structure containing wavenumber information and coherence values.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `N`: Number of wavenumber dimensions
- `A`: Type of wavenumber argument
- `T`: Type of coherence values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Fields
- `wavenumber`: Wavenumber grid or values
- `coherence`: Coherence estimates
- `processinformation`: Information about the processes
- `estimationinformation`: Information about the estimation procedure
"""
struct Coherence{E, D, N, A, T, IP, IE} <: AbstractEstimate{E, D, N}
    wavenumber::A
    coherence::T
    processinformation::IP
    estimationinformation::IE
    function Coherence{E}(
            wavenumber::NTuple{N}, coherence::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E <: EstimateTrait, D, N, T, IE}
        checkinputs(wavenumber, coherence, processinfo)
        IP = typeof(processinfo)
        A = typeof(wavenumber)
        return new{E, D, N, A, T, IP, IE}(
            wavenumber, coherence, processinfo, estimationinfo)
    end
    function Coherence{E}( # for inputs that are rotational
            wavenumber::A, coherence::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E <: EstimateTrait, D, A, T, IE}
        checkinputs(wavenumber, coherence, processinfo)
        IP = typeof(processinfo)
        return new{E, D, 1, A, T, IP, IE}(
            wavenumber, coherence, processinfo, estimationinfo)
    end
end
const RotationalCoherence{E, D, S <: Coherence} = RotationalEstimate{E, D, S}
const NormalOrRotationalCoherence{E} = Union{Coherence{E}, RotationalCoherence{E}}

## required interface
# preference given to direct computation from marginal spectra in partial case
# computation from partial spectra will give different results in general due to debiasing
# This debiasing is appropriate for linear transformations of the partial spectra, which is
# why is is not the default for coherence computation
computed_from(::Type{<:Coherence{MarginalTrait, D}}) where {D} = Spectra{MarginalTrait, D}
function computed_from(::Type{<:Coherence{PartialTrait, D}}) where {D}
    (Spectra{MarginalTrait, D}, Spectra{PartialTrait, D})
end
function computed_from(::Type{<:RotationalCoherence{MarginalTrait, D}}) where {D}
    RotationalSpectra{MarginalTrait, D}
end
function computed_from(::Type{<:RotationalCoherence{PartialTrait, D}}) where {D}
    (RotationalSpectra{MarginalTrait, D}, RotationalSpectra{PartialTrait, D})
end

function allocate_estimate_memory(::Type{<:NormalOrRotationalCoherence},
        ::Type{<:NormalOrRotationalSpectra}, relevant_memory; kwargs...)
    return relevant_memory, nothing
end

function extract_relevant_memory(::Type{<:NormalOrRotationalCoherence},
        est::NormalOrRotationalSpectra)
    return get_estimates(est)
end
function extract_relevant_memory(::Type{<:NormalOrRotationalCoherence},
        est::EstimateMemory{<:NormalOrRotationalSpectra})
    return est.output_memory
end

function validate_core_parameters(
        ::Type{<:NormalOrRotationalCoherence}; kwargs...)
    # no additional parameters to validate
    return nothing
end

function resolve_missing_parameters(
        ::Type{<:NormalOrRotationalCoherence}, arg; kwargs...)
    return kwargs
end

function validate_memory_compatibility(
        ::Type{<:NormalOrRotationalCoherence}, mem, source; kwargs...)
    # no additional compatibility checks needed
    @argcheck size(mem.output_memory) == size(get_estimates(source))
    @argcheck eltype(mem.output_memory) == eltype(get_estimates(source))
end

function compute_estimate!(::Type{<:NormalOrRotationalCoherence{E}},
        mem, source::Spectra{E}; kwargs...) where {E}
    if !is_same_process_sets(source)
        throw(ArgumentError(
            "Coherence computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(source))."
        ))
    end

    return _compute_coherence_estimate!(_coherence_noalloc!, mem.output_memory, source, E)
end

function compute_estimate!(::Type{<:NormalOrRotationalCoherence{PartialTrait}},
        mem, source::Spectra{MarginalTrait}; kwargs...)
    if !is_same_process_sets(source)
        throw(ArgumentError(
            "Coherence computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(source))."
        ))
    end

    return _compute_coherence_estimate!(
        _partial_coherence_noalloc!, mem.output_memory, source, PartialTrait)
end

get_evaluation_points(est::Coherence) = est.wavenumber

get_estimates(est::Coherence) = est.coherence

## internals

function coherence(x::SMatrix)
    d = diagm(sqrt.(inv.(diag(x))))
    return d * x * d
end
function coherence(x::AbstractMatrix)
    y = deepcopy(x)
    return coherence!(y)
end
function coherence!(x::AbstractMatrix)
    for i in axes(x, 1), j in axes(x, 2)
        if i == j
            continue
        end
        x[i, j] /= sqrt(x[i, i] * x[j, j])
    end
    for i in axes(x, 1)
        x[i, i] = one(eltype(x))
    end
    return x
end

coherence(x::Number) = one(typeof(x))

_coherence_noalloc!(x::Union{Number, SMatrix}) = coherence(x)
_coherence_noalloc!(x::AbstractMatrix) = coherence!(x)

function partial_coherence(x::SMatrix)
    return -coherence(inv(x)) + 2I # add 2I to set diagonals to 1
end
function partial_coherence(x::AbstractMatrix)
    y = deepcopy(x)
    return partial_coherence!(y)
end
function partial_coherence!(x::AbstractMatrix)
    x = -coherence!(LinearAlgebra.inv!(cholesky!(x)))
    for i in axes(x, 1)
        x[i, i] = one(eltype(x))
    end
    return x
end

partial_coherence(x::Number) = one(typeof(x))

_partial_coherence_noalloc!(x::Union{Number, SMatrix}) = partial_coherence(x)
_partial_coherence_noalloc!(x::AbstractMatrix) = partial_coherence!(x)

"""
    _compute_coherence_estimate!(spectrum, transform_func!, trait_type)

Internal helper to compute coherence estimates with proper error handling.
"""
function _compute_coherence_estimate!(transform_func!, value, spectrum, trait_type)
    trait = process_trait(spectrum)
    est = get_estimates(spectrum)
    process_info = get_process_information(spectrum)
    estimation_info = get_estimation_information(spectrum)
    value = apply_transform!(transform_func!, value, est, trait)
    return Coherence{trait_type}(
        get_evaluation_points(spectrum), value, process_info, estimation_info)
end
