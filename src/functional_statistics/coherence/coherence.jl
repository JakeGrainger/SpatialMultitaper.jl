"""
    coherence(spectrum::NormalOrRotationalSpectra{E}) where {E} -> Coherence{E}
    coherence(data, region; kwargs...) -> Coherence
    coherence(data::SpatialData; kwargs...) -> Coherence

Compute coherence from a spectral estimate or directly from spatial data.

The coherence function computes the normalized cross-spectral density, providing a measure
of linear dependence between processes as a function of wavenumber. For a cross-spectral
matrix S, the coherence γᵢⱼ is computed as γᵢⱼ = Sᵢⱼ / √(Sᵢᵢ * Sⱼⱼ).

# Arguments
- `spectrum::NormalOrRotationalSpectra{E}`: A `Spectra` or `RotationalSpectra` estimate
    with multiple processes
- `data`: Spatial data for direct coherence computation
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
- `Coherence{E}`: A coherence estimate object containing:
  - `wavenumber`: Wavenumber grid matching the input spectrum
  - `coherence`: Coherence estimates with the same spatial dimensions as input
  - `processinformation`: Information about the analyzed processes
  - `estimationinformation`: Details about the estimation procedure

# Throws
- `ArgumentError`: If the spectrum does not have equal input and output process sets
    (i.e., the spectral matrix is not square)

# Notes
- Coherence values are complex numbers with magnitude ≤ 1
- Diagonal elements are always 1 (perfect self-coherence)
- For single-process data, returns scalar coherence of 1
- Use [`magnitude_coherence`](@ref) or [`magnitude_squared_coherence`](@ref) for
    magnitude-only results

# Examples
```julia
# Compute coherence from existing spectral estimate
spec = spectra(data; kmax = 0.5, nw = 3)
coh = coherence(spec)

# Direct computation from data and region
coh = coherence(data, region; nk = (32, 32), kmax = (0.5, 0.5), nw = 4)

# Direct computation from SpatialData object
spatial_data_obj = spatial_data(data, region)
coh = coherence(spatial_data_obj; kmax = 0.3, tapers = my_tapers)

# Access coherence between processes 1 and 2
cross_coherence = coh[1, 2]
```
"""
coherence

"""
    partial_coherence(spectrum) -> Coherence{PartialTrait}
    partial_coherence(data, region; kwargs...) -> Coherence{PartialTrait}
    partial_coherence(data::SpatialData; kwargs...) -> Coherence{PartialTrait}

Compute partial coherence from a spectral estimate or directly from spatial data.

See [`coherence`](@ref) for details on arguments, keywords, and return types.
"""
partial_coherence

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
get_evaluation_points(est::Coherence) = est.wavenumber
get_estimates(est::Coherence) = est.coherence

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

function coherence(spectrum::NormalOrRotationalSpectra{E})::Coherence{E} where {E}
    mem = deepcopy(spectrum)
    return coherence!(mem)
end
function coherence!(spectrum::NormalOrRotationalSpectra{E})::Coherence{E} where {E}
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Coherence computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(spectrum))."
        ))
    end

    return _compute_coherence_estimate!(spectrum, _coherence_noalloc!, E)
end
_coherence_noalloc!(x::Union{Number, SMatrix}) = coherence(x)
_coherence_noalloc!(x::AbstractMatrix) = coherence!(x)

function coherence(data, region; kwargs...)::Coherence
    return coherence(spatial_data(data, region); kwargs...)
end

function coherence(data::SpatialData; kwargs...)::Coherence
    return coherence(spectra(data; kwargs...))
end

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

function partial_coherence(spectrum::NormalOrRotationalSpectra)::Coherence{PartialTrait}
    mem = deepcopy(spectrum)
    return partial_coherence!(mem)
end

function partial_coherence!(spectrum::NormalOrRotationalSpectra{PartialTrait})::Coherence{PartialTrait}
    return coherence!(spectrum) # partial coherence is just coherence of the partial spectra
end

function partial_coherence!(spectrum::NormalOrRotationalSpectra{MarginalTrait})
    if !is_same_process_sets(spectrum)
        throw(ArgumentError(
            "Partial coherence computation requires equal input and output process sets. " *
            "Got processes with dimensions $(size(spectrum))."
        ))
    end

    return _compute_coherence_estimate!(spectrum, _partial_coherence_noalloc!, PartialTrait)
end
_partial_coherence_noalloc!(x::Union{Number, SMatrix}) = partial_coherence(x)
_partial_coherence_noalloc!(x::AbstractMatrix) = partial_coherence!(x)

function partial_coherence(data, region; kwargs...)::Coherence{PartialTrait}
    return partial_coherence(spatial_data(data, region); kwargs...)
end

function partial_coherence(data::SpatialData; kwargs...)::Coherence{PartialTrait}
    return partial_coherence(spectra(data; kwargs...))
end

"""
    magnitude_coherence(spectrum::Spectra)

Compute magnitude coherence |γᵢⱼ| from a spectral estimate.
"""
magnitude_coherence(spectrum::Spectra) = abs(coherence(spectrum))

"""
    magnitude_coherence(coh::Coherence)

Compute magnitude of existing coherence estimate.
"""
magnitude_coherence(coh::Coherence) = abs(coh)

"""
    magnitude_coherence(args...; kwargs...)

Compute magnitude coherence directly from data arguments.
"""
function magnitude_coherence(args...; kwargs...)
    return magnitude_coherence(spectra(args...; kwargs...))
end

"""
    partial_magnitude_coherence(spectrum::Spectra{MarginalTrait})

Compute partial magnitude coherence from marginal spectral estimates.
"""
function partial_magnitude_coherence(spectrum::Spectra{MarginalTrait})
    return magnitude_coherence(spectrum)
end

function partial_magnitude_coherence(spectrum::Spectra{PartialTrait})
    return magnitude_coherence(partial_spectra(spectrum))
end

function partial_magnitude_coherence(coh::Coherence{MarginalTrait})
    throw(_partial_from_marginal_error("partial magnitude coherence", typeof(coh)))
end

partial_magnitude_coherence(coh::Coherence{PartialTrait}) = magnitude_coherence(coh)

function partial_magnitude_coherence(args...; kwargs...)
    return partial_magnitude_coherence(spectra(args...; kwargs...))
end

"""
    magnitude_squared_coherence(spectrum::Spectra)

Compute magnitude-squared coherence |γᵢⱼ|² from a spectral estimate.

Also known as the coherency or coherence function.
"""
magnitude_squared_coherence(spectrum::Spectra) = abs2(coherence(spectrum))

magnitude_squared_coherence(coh::Coherence) = abs2(coh)

function magnitude_squared_coherence(args...; kwargs...)
    return magnitude_squared_coherence(spectra(args...; kwargs...))
end

function partial_magnitude_squared_coherence(spectrum::Spectra{MarginalTrait})
    return magnitude_squared_coherence(spectrum)
end

function partial_magnitude_squared_coherence(spectrum::Spectra{PartialTrait})
    return magnitude_squared_coherence(partial_spectra(spectrum))
end

function partial_magnitude_squared_coherence(coh::Coherence{MarginalTrait})
    throw(_partial_from_marginal_error("partial magnitude squared coherence", typeof(coh)))
end

function partial_magnitude_squared_coherence(coh::Coherence{PartialTrait})
    return magnitude_squared_coherence(coh)
end

function partial_magnitude_squared_coherence(args...; kwargs...)
    return partial_magnitude_squared_coherence(spectra(args...; kwargs...))
end

"""
    phase(spectrum::Union{Spectra, Coherence})

Compute phase of coherence or cross-spectral estimates.

Returns the phase angle in radians.
"""
phase(spectrum::Union{Spectra, Coherence}) = angle(spectrum)

phase(args...; kwargs...) = phase(spectra(args...; kwargs...))

function partial_phase(spectrum::Union{Spectra{PartialTrait}, Coherence{PartialTrait}})
    return phase(spectrum)
end

partial_phase(spectrum::Spectra{MarginalTrait}) = phase(partial_spectra(spectrum))

function partial_phase(coh::Coherence{MarginalTrait})
    throw(_partial_from_marginal_error("partial phase", typeof(coh)))
end

partial_phase(args...; kwargs...) = partial_phase(spectra(args...; kwargs...))

# Helper functions

"""
    _compute_coherence_estimate!(spectrum, transform_func, trait_type)

Internal helper to compute coherence estimates with proper error handling.
"""
function _compute_coherence_estimate!(spectrum, transform_func, trait_type)
    trait = process_trait(spectrum)
    est = get_estimates(spectrum)
    process_info = get_process_information(spectrum)
    estimation_info = get_estimation_information(spectrum)
    transformed = apply_transform!(transform_func, est, trait)
    return Coherence{trait_type}(
        get_evaluation_points(spectrum), transformed, process_info, estimation_info)
end

"""
    _partial_from_marginal_error(operation_name, estimate_type)

Generate a descriptive error for unsupported partial operations on marginal estimates.
"""
function _partial_from_marginal_error(operation_name::String, estimate_type::Type)
    return ArgumentError(
        "Cannot compute $(operation_name) from a marginal coherence estimate. " *
        "Compute from spectral estimates or use partial spectral estimates instead. " *
        "Got estimate of type: $(estimate_type)"
    )
end
