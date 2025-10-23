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
    return (mem = mem, power = power)
end

extract_relevant_memory(::Type{<:Spectra{MarginalTrait}}, source) = nothing

function validate_core_parameters(::Type{<:Spectra{MarginalTrait}}, kwargs...)
end

# function validate_memory_compatibility end

function resolve_missing_parameters(
        ::Type{<:Spectra{MarginalTrait}}, data::SpatialData, kwargs...)
end
# function compute_estimate! end

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

function validate_dk end
function validate_nk end
function validate_kmax end
function resolve_dk_nk_kmax end
