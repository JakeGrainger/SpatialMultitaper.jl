struct Spectra{E, D, A, T, IP, IE} <: AnisotropicEstimate{E, D}
    wavenumber::NTuple{D, A}
    power::T
    processinformation::IP
    estimationinformation::IE
    function Spectra{E}(wavenumber::NTuple{D, A}, power::T, processinfo::IP,
            estimationinfo::IE) where {E <: EstimateTrait, D, A, T, IP, IE}
        checkinputs(wavenumber, power, processinfo)
        return new{E, D, A, T, IP, IE}(wavenumber, power, processinfo, estimationinfo)
    end
end

const RotationalSpectra{E, D, S <: Spectra} = RotationalEstimate{E, D, S}
const NormalOrRotationalSpectra{E} = Union{Spectra{E}, RotationalSpectra{E}}
get_short_estimate_name(::Type{<:Spectra}) = "f"

## interface
computed_from(::Type{<:Spectra{MarginalTrait}}) = SpatialData

function allocate_estimate_memory(::Type{<:Spectra{MarginalTrait}}, ::Type{T},
        ::Nothing; nk, kmax, tapers, kwargs...) where {T <: SpatialData}
    mem = preallocate_tapered_dft(T, tapers, nk, kmax)
    power = _preallocate_spectral_matrix(T, nk)
    return power, mem
end

extract_relevant_memory(::Type{<:Spectra{MarginalTrait}}, source) = nothing

function validate_core_parameters(::Type{<:Spectra{MarginalTrait}}, kwargs...)
    spectral_kwargs = (:nk, :kmax, :dk, :tapers, :nw, :mean_method)
    for (k, v) in kwargs
        if k in spectral_kwargs
            validate_spectral_kwarg(k, v)
        end
    end
    return nothing
end

function validate_spectral_kwarg(k::Symbol, v)
    if k == :nk
        validate_nk(v)
    elseif k == :kmax
        validate_kmax(v)
    elseif k == :dk
        validate_dk(v)
    elseif k == :tapers
        validate_tapers(v)
    elseif k == :nw
        validate_nw(v)
    elseif k == :mean_method
        validate_mean_estimation_method(v)
    else
        error("Unknown spectral kwarg: $k")
    end
end

function resolve_missing_parameters(
        ::Type{<:Spectra{MarginalTrait}}, data::SpatialData, kwargs...)
    stage_1 = resolve_dk_nk_kmax(data; kwargs...)
    stage_2 = resolve_tapers(data; stage_1...)
    return stage_2
end

function validate_memory_compatibility(
        ::Type{<:Spectra{MarginalTrait}}, mem, source; kwargs...)
    validate_dft_memory(mem.internal_memory, source; kwargs...)
    validate_spectral_memory(mem.output_memory, source; kwargs...)
end

function compute_estimate!(::Type{<:Spectra{MarginalTrait}}, mem, source::SpatialData;
        tapers, nk, kmax, mean_method, kwargs...)
    dft_mem = mem.internal_memory
    power = mem.output_memory
    J_n = tapered_dft!(dft_mem, data, tapers, nk, kmax, mean_method)
    _dft_to_spectral_matrix!(power, J_n, process_trait(data))

    process_info = ProcessInformation(data; mean_method = mean_method)
    estimation_info = EstimationInformation(length(tapers))

    return Spectra{MarginalTrait}(wavenumber, power, process_info, estimation_info)
end

get_evaluation_points(est::Spectra) = est.wavenumber

get_estimates(est::Spectra) = est.power

## Internals

### preallocation

function _preallocate_spectral_matrix(::Type{<:SingleProcessData}, nk)
    return zeros(Float64, nk) # TODO: generalize type
end

function _preallocate_spectral_matrix(::Type{<:MultipleSpatialDataTuple{P}}, nk) where {P}
    return zeros(SMatrix{P, P, ComplexF64, P * P}, nk) # TODO: generalize type
end

function _preallocate_spectral_matrix(::Type{<:MultipleSpatialDataVec}, nk)
    return zeros(ComplexF64, (ncol(data), ncol(data), nk...)) # TODO: generalize type
end

### keywords
include("parameter_validation.jl")

### computation
include("computation.jl")

### memory validation

function validate_dft_memory(mem::AbstractArray, source; kwargs...)
    # TODO: implement this
end

function validate_spectral_memory(
        mem::AbstractArray, source::SingleProcessData; nk, kwargs...)
    @argcheck size(mem) == nk
    @argcheck eltype(mem) == ComplexF64 # TODO: generalize type
end
function validate_spectral_memory(
        mem::AbstractArray, source::MultipleSpatialDataTuple{P}; nk, kwargs...) where {P}
    @argcheck size(mem) == nk
    @argcheck eltype(mem) <: SMatrix{P, P, ComplexF64} # TODO: generalize type
end
function validate_spectral_memory(
        mem::AbstractArray, source::MultipleSpatialDataVec; nk, kwargs...)
    nprocesses = ncol(source)
    @argcheck size(mem) == (nprocesses, nprocesses, nk)
    @argcheck eltype(mem) == ComplexF64 # TODO: generalize type
end
