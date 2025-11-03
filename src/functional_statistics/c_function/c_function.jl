"""
    CFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}

Spatial C function estimate derived from spectral estimates.

The C function represents the reduced covariance measure of a stationary process evaluated
on punctured balls as a function of distance. At radius `r`, it computes
C(r) = C̆({x : 0 < ||x|| ≤ r}), where C̆ is the reduced covariance measure. This provides
a spatial domain representation of correlation structure complementary to spectral estimates.

# Type Parameters
- `E <: EstimateTrait`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension of the underlying process
- `A`: Type of the radii array
- `T`: Type of the C function values (typically `Vector{Float64}` or similar)
- `IP`: Type of process information structure
- `IE`: Type of estimation information structure

# Fields
- `radii::A`: Distance values at which the C function is evaluated
- `value::T`: C function values corresponding to each radius
- `processinformation::IP`: Information about the analyzed processes
- `estimationinformation::IE`: Details about the estimation procedure

# Mathematical Background
For a stationary process with power spectral density f(k), the C function is computed as:

C(r) = ∫ f(k) W(r,k) dk

where W(r,k) is the spatial weighting function determined by the dimension D and represents
the Fourier transform of the characteristic function of a ball of radius r.

# Notes
- C function values are always real-valued due to spatial integration
- The function measures cumulative spatial correlation within distance r
- Supports both marginal and partial variants via the trait system
- Integrates over punctured balls (excluding the origin) to avoid singularities

See also: [`c_function`](@ref), [`partial_c_function`](@ref)
"""
struct CFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function CFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        return new{E, D, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end
get_short_base_estimate_name(::Type{<:CFunction}) = "C"
get_base_estimate_name(::Type{<:CFunction}) = "C function"

## required interface

function computed_from(::Type{<:CFunction{E, D}}) where {E, D}
    return (Spectra{E, D}, RotationalSpectra{E, D})
end

function select_source_type(
        T::Type{<:CFunction}, arg; c_algorithm = CFunctionAutoSelect(), kwargs...)
    _select_c_source_type(T, c_algorithm, arg; kwargs...)
end

abstract type CFunctionSelection end
struct CFunctionFromSpectra <: CFunctionSelection end
struct CFunctionFromRotationalSpectra <: CFunctionSelection end
struct CFunctionAutoSelect <: CFunctionSelection end
function _select_c_source_type(
        ::Type{<:CFunction{E, D}}, ::CFunctionFromSpectra, arg; kwargs...) where {E, D}
    return Spectra{E, D}
end
function _select_c_source_type(::Type{<:CFunction{E, D}}, ::CFunctionFromRotationalSpectra,
        arg; kwargs...) where {E, D}
    return RotationalSpectra{E, D}
end
function _select_c_source_type(
        T::Type{<:CFunction{E, D}}, ::CFunctionAutoSelect, arg; kwargs...) where {E, D}
    return Spectra{E, D}
end
function _select_c_source_type(T::Type{<:CFunction{E, D}}, ::CFunctionAutoSelect,
        arg::RotationalSpectra{E, D}; kwargs...) where {
        E, D}
    return RotationalSpectra{E, D}
end
process_c_algorithm(s::String) = process_c_algorithm(Symbol(s))
process_c_algorithm(s::Symbol) = process_c_algorithm(Val{s}())
process_c_algorithm(::Nothing) = CFunctionAutoSelect()
process_c_algorithm(::Val{:auto}) = CFunctionAutoSelect()
process_c_algorithm(::Val{:from_spectra}) = CFunctionFromSpectra()
process_c_algorithm(::Val{:from_rotational_spectra}) = CFunctionFromRotationalSpectra()
function process_c_algorithm(::Val{T}) where {T}
    throw(ArgumentError("Unknown C function algorithm selection: $T"))
end

function _extract_wavenumber_from_c_mem(::Type{<:Spectra}, ::Nothing; nk, kmax, kwargs...)
    _choose_wavenumbers_1d.(nk, kmax)
end
function _extract_wavenumber_from_c_mem(
        ::Type{<:RotationalSpectra}, ::Nothing; rotation_radii, kwargs...)
    rotation_radii
end
function _extract_wavenumber_from_c_mem(
        ::Type{<:NormalOrRotationalSpectra}, wavenumber; kwargs...)
    wavenumber
end

function allocate_estimate_memory(::Type{<:CFunction}, ::Type{S}, relevant_memory;
        kwargs...) where {S <: NormalOrRotationalSpectra}
    mem = relevant_memory[1:2]
    wavenumber = _extract_wavenumber_from_c_mem(S, relevant_memory[3]; kwargs...)
    spatial_output = preallocate_c_output(S, mem...; kwargs...)
    weights = precompute_c_weights(S, mem[1], wavenumber; kwargs...)
    return spatial_output, weights
end
function preallocate_c_output(::Type{<:NormalOrRotationalSpectra}, mem...; kwargs...)
    return preallocate_radial_output(mem...; kwargs...)
end

function extract_relevant_memory(::Type{<:CFunction}, est::NormalOrRotationalSpectra)
    return deepcopy(get_estimates(est)), process_trait(est), get_evaluation_points(est)
end
function extract_relevant_memory(
        ::Type{<:CFunction}, mem::EstimateMemory{<:NormalOrRotationalSpectra})
    return mem.output_memory, process_trait(mem), nothing
end

function preprocess_algorithm_kwargs(::Type{<:CFunction}; kwargs...)
    if :c_algorithm in keys(kwargs)
        c_algorithm = process_c_algorithm(kwargs[:c_algorithm])
        other_kwargs = filter(kv -> kv[1] != :c_algorithm, kwargs)
        return (c_algorithm = c_algorithm, other_kwargs...)
    else
        return kwargs
    end
end

function validate_core_parameters(::Type{<:CFunction}; kwargs...)
    if :radii in keys(kwargs)
        radii = kwargs[:radii]
        validate_radii(radii)
    end
    return nothing
end

function resolve_missing_parameters(::Type{<:CFunction}, arg; kwargs...)
    radii = get(kwargs, :radii, nothing)
    c_algorithm = get(kwargs, :c_algorithm, CFunctionAutoSelect()) # Already processed, use default if missing
    other_kwargs = filter(kv -> kv[1] ∉ (:radii, :c_algorithm), kwargs)
    return (radii = process_radii(radii, arg),
        c_algorithm = c_algorithm, other_kwargs...)
end

function validate_memory_compatibility(
        ::Type{<:CFunction}, mem, arg::NormalOrRotationalSpectra; radii, kwargs...)
    validate_c_internal(mem.internal_memory, get_estimates(arg), process_trait(arg))
    validate_radial_memory(mem.output_memory, process_trait(arg), radii)
    return nothing
end

function compute_estimate!(
        ::Type{<:CFunction{E}}, mem, source::NormalOrRotationalSpectra; radii, kwargs...) where {E}
    out = mem.output_memory
    store = mem.internal_memory
    value = _sdf2C!(out, store, source, radii)

    process_info = get_process_information(source)
    estimation_info = get_estimation_information(source)
    return CFunction{E}(radii, value, process_info, estimation_info)
end

get_evaluation_points(f::CFunction) = f.radii

get_estimates(f::CFunction) = f.value

## Internal

### Validation

function validate_c_internal(weights, power, ::SingleProcessTrait)
    @argcheck size(weights)[1:(end - 1)] == size(power)
    return nothing
end

function validate_c_internal(weights, power, ::MultipleTupleTrait)
    @argcheck size(weights)[1:(end - 1)] == size(power)
    return nothing
end

function validate_c_internal(weights, power, ::MultipleVectorTrait)
    @argcheck size(weights)[1:(end - 1)] == size(power)[3:end]
    return nothing
end

### Computation

function _sdf2C!(out, store, spectrum::NormalOrRotationalSpectra, radii)
    power = get_estimates(spectrum)
    zero_atom = get_process_information(spectrum).atoms
    trait = process_trait(spectrum)
    return _sdf2C_internal!(out, store, power, trait, zero_atom, radii)
end

function _sdf2C_internal!(
        out, store, power, ::Union{SingleProcessTrait, MultipleTupleTrait},
        zero_atom, radii::AbstractVector)
    for i in eachindex(out)
        out[i] = _sdf2C_sum(selectdim(store, ndims(store), i), power, zero_atom)
    end
    return out
end

function _sdf2C_internal!(
        out, store, power::AbstractArray{<:Number, N}, ::MultipleVectorTrait,
        zero_atom, radii::AbstractVector) where {N}
    for i in axes(out, ndims(out))
        store_slice = selectdim(store, ndims(store), i)
        for idx in CartesianIndices(size(out)[1:2])
            power_slice = view(power, idx, ntuple(Returns(:), N - 2)...)
            zero_atom_slice = _slice_zero_atom(zero_atom, idx)
            out[idx, i] = _sdf2C_sum(store_slice, power_slice, zero_atom_slice)
        end
    end
    return out
end

_slice_zero_atom(zero_atom, idx) = zero_atom[idx]
_slice_zero_atom(::Nothing, _) = nothing

# Helper functions for C function computation
function _sdf2C_sum(store, power, zero_atom)
    c = sum((p - zero_atom) * s for (p, s) in zip(power, store))
    return real(c)
end
function _sdf2C_sum(store, power, ::Nothing)
    c = sum(p * s for (p, s) in zip(power, store))
    return real(c)
end

function precompute_c_weights(::Type{<:Spectra{E, D}}, power, wavenumber;
        radii::AbstractVector, kwargs...) where {E, D}
    wavenumber_size = length.(wavenumber)
    wavenumber_spacing = prod(step, wavenumber)
    T = promote_type(float(eltype(radii)), float(typeof(wavenumber_spacing)))
    store = zeros(T, wavenumber_size..., length(radii))
    for (i, radius) in enumerate(radii)
        for (idx, k) in zip(
            CartesianIndices(wavenumber_size), Iterators.product(wavenumber...))
            store[idx, i] = _anisotropic_c_weight(radius, k, Val{D}()) * wavenumber_spacing
        end
    end
    return store
end

function precompute_c_weights(::Type{<:RotationalSpectra{E, D}}, power,
        wavenumber; radii::AbstractVector, kwargs...) where {E, D}
    spacing = step(wavenumber)
    T = promote_type(float(eltype(radii)), float(typeof(spacing)))
    store = zeros(T, length(wavenumber), length(radii))
    for (i, radius) in enumerate(radii)
        for (idx, k) in enumerate(wavenumber)
            store[idx, i] = _isotropic_c_weight(radius, k, spacing, Val{D}())
        end
    end
    return store
end

# Spatial weighting functions

"""
    _anisotropic_c_weight(r, u, ::Val{D})

Compute spatial weighting function for D-dimensional C function computation.

These functions represent the Fourier transform of the characteristic function of a
D-dimensional ball of radius r, evaluated at wavenumber u. They are the fundamental
building blocks for transforming spectral estimates to spatial C functions.

# Mathematical Details
- **1D**: W(r,|u|) = 2r sinc(2r|u|)
- **2D**: W(r,|u|) = (r/|u|) J₁(2πr|u|) for |u| > 0, πr² for u = 0
- **3D**: W(r,|u|) = (r/|u|)^(3/2) J₃/₂(2πr|u|) for |u| > 0, (4π/3)r³ for u = 0
- **General D**: Uses Bessel functions J_{D/2} with appropriate normalization

where J_ν denotes the Bessel function of the first kind of order ν.
"""
function _anisotropic_c_weight(r, u, ::Val{1})
    x = norm(u)
    return 2r * sinc(2r * x) # sinc(x) = sin(πx)/(πx)
end

function _anisotropic_c_weight(r, u, ::Val{2})
    x = norm(u)
    return (x < 1e-10) ? (π * r^2) : ((r / x) * besselj1(2π * r * x))
end

function _anisotropic_c_weight(r, u, ::Val{3})
    x = norm(u)
    return (x < 1e-10) ? (4π / 3 * r^3) : (sqrt(r / x)^3 * besselj(3 / 2, 2π * r * x))
end

function _anisotropic_c_weight(r, u, ::Val{D}) where {D}
    x = norm(u)
    return _anisotropic_c_weight_generic(r, x, Val{D}())
end

function _anisotropic_c_weight_generic(r, x, ::Val{D}) where {D}
    if x < 1e-10
        # Handle singularity at origin using ball measure
        return unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), r))
    else
        return (r / x)^(D / 2) * besselj(D / 2, 2π * r * x)
    end
end

"""
    _isotropic_c_weight(r, k, spacing, ::Val{D})

Compute weighting function for isotropic C functions from rotational spectra.

For rotationally averaged spectral estimates, this function computes the appropriate
weights that account for the radial integration already performed in the spectral
domain. The spacing parameter accounts for the discrete wavenumber grid integration.

# Arguments
- `r`: Spatial radius for C function evaluation
- `k`: Radial wavenumber coordinate
- `spacing`: Wavenumber spacing for integration
- `::Val{D}`: Spatial dimension

# Mathematical Details
The weights integrate the anisotropic weights over annular regions in wavenumber space,
accounting for the discrete grid spacing used in the rotational averaging process.
"""
function _isotropic_c_weight(r, k, spacing, ::Val{1})
    half_spacing = spacing / 2
    # note sinint is for the unnormalised sinc integral
    # note 2 times because we integrate from k - half_spacing to k + half_spacing and from
    # -k - half_spacing to -k + half_spacing in this case
    return 2 / π *
           (sinint(2π * r * (k + half_spacing)) - sinint(2π * r * (k - half_spacing)))
end
function _isotropic_c_weight(r, k, spacing, ::Val{2})
    half_spacing = spacing / 2
    return besselj0(2π * r * (k - half_spacing)) - besselj0(2π * r * (k + half_spacing))
end

function _isotropic_c_weight(r, k, spacing, ::Val{3})
    half_spacing = spacing / 2
    return 2 / π * (sinint(2π * r * (k + half_spacing)) - sin(2π * r * (k + half_spacing)) -
            sinint(2π * r * (k - half_spacing)) + sin(2π * r * (k - half_spacing)))
end

function _isotropic_c_weight(r, k, spacing, ::Val{D}) where {D}
    _isotropic_c_weight_generic(r, k, spacing, Val{D}())
end

function _isotropic_c_weight_generic(r, k, spacing, ::Val{D}) where {D}
    half_spacing = spacing / 2
    _isotropic_c_weight_generic_int(r, k + half_spacing, D) -
    _isotropic_c_weight_generic_int(r, k - half_spacing, D)
end

function _isotropic_c_weight_generic_int(r, x, d)
    u = pi * r * x
    return u^d / (gamma(d / 2 + 1)^2) * pFq((d / 2,), (d / 2 + 1, d / 2 + 1), -u^2)
end
