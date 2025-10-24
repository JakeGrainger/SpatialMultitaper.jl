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

function computed_from(::Type{CFunction{E}}) where {E}
    return (Spectra{E}, RotationalSpectra{E})
end

function allocate_estimate_memory(
        ::Type{<:CFunction}, ::Type{<:Spectra}, previous_memory; radii, kwargs...)
    spatial_output = preallocate_spatial_output(spectrum, radii)
    weights = precompute_c_weights(spectrum, radii)
    return (spatial_output = spatial_output, weights = weights)
end

function extract_relevant_memory(::Type{CFunction}, est::Spectra)
    return get_estimates(est)
end
function extract_relevant_memory(::Type{CFunction}, mem::EstimateMemory{<:Spectra})
    return mem.power
end

function validate_core_parameters(::Type{<:CFunction}, radii, kwargs...)
    @argcheck all(radii .>= 0)
    side_length = sides(boundingbox(getregion(data)))
    @argcheck all(radii .<= Meshes.ustrip(minimum(side_length)))
    return nothing
end
# when no radii are provided, defaults will be set in apply_parameter_defaults
validate_core_parameters(::Type{<:CFunction}, kwargs...) = nothing

function resolve_missing_parameters(::Type{<:CFunction}, data::SpatialData; kwargs...)
    return (radii = process_radii(radii, data), kwargs...)
end

process_radii(radii, ::SpatialData) = radii
function process_radii(::Nothing, data::SpatialData)
    region = getregion(data)
    short_side = Meshes.ustrip(minimum(sides(boundingbox(region))))
    return range(0, short_side / 3, length = 100)
end

function validate_memory_compatibility(
        ::Type{<:CFunction}, mem, arg::Spectra; radii, kwargs...)
    validate_c_internal(mem.weights, get_estimates(arg), get_trait(arg))
    validate_c_output(mem.spatial_output, radii, get_trait(arg))
    return nothing
end

function compute_estimate! end

get_evaluation_points(f::CFunction) = f.radii

get_estimates(f::CFunction) = f.value

## Internal

### Validation

function validate_c_internal(weights, spectrum, ::SingleProcessTrait)
    @argcheck size(weights) == size(spectrum)
    return nothing
end

function validate_c_internal(weights, spectrum, ::MultipleSpatialDataTuple)
    @argcheck size(weights) == size(spectrum)
    return nothing
end

function validate_c_internal(weights, spectrum, ::MultipleSpatialDataVec)
    @argcheck size(weights) == size(spectrum)[3:end]
    return nothing
end

function validate_c_output(output::AbstractVector, radii, ::SingleProcessTrait)
    @argcheck length(output) == length(radii)
    @argcheck eltype(output) <: Real
    return nothing
end

function validate_c_output(output::AbstractVector, radii, ::MultipleSpatialDataTuple)
    @argcheck length(output) == length(radii)
    @argcheck eltype(output) <: SMatrix
    return nothing
end

function validate_c_output(output::AbstractArray, radii, ::MultipleSpatialDataVec)
    @argcheck ndims(output) == 3
    @argcheck size(output, 3) == length(radii)
    @argcheck eltype(output) <: Real
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

function precompute_c_weights(spectra::Spectra{E, D}, radii::AbstractVector) where {E, D}
    wavenumber = get_evaluation_points(spectra)
    power_size = size(get_estimates(spectra))
    wavenumber_spacing = prod(step, wavenumber)
    T = promote_type(float(eltype(radii)), float(typeof(wavenumber_spacing)))
    store = zeros(T, power_size..., length(radii))
    for (i, radius) in enumerate(radii)
        for (idx, k) in zip(CartesianIndices(power_size), Iterators.product(wavenumber...))
            store[idx, i] = _anisotropic_c_weight(radius, k, Val{D}()) * wavenumber_spacing
        end
    end
    return store
end

function precompute_c_weights(
        spectra::RotationalSpectra{E, D}, radii::AbstractVector) where {E, D}
    wavenumber = get_evaluation_points(spectra)
    power_size = size(get_estimates(spectra))
    spacing = step(wavenumber)
    T = promote_type(float(eltype(radii)), float(typeof(spacing)))
    store = zeros(T, power_size..., length(radii))
    for (i, radius) in enumerate(radii)
        for (idx, k) in zip(CartesianIndices(power_size), wavenumber)
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
