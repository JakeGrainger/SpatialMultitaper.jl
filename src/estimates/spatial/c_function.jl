"""
    CFunction{E, D, A, T, IP, IE} <: IsotropicEstimate{E, D}

Spatial C function estimate derived from spectral estimates.

The C function is the reduced covariance measure of a process evaluated on a punctured ball,
of a function of distance. So at radius `r` it is C(r) = C̆({x : 0 < ||x|| ≤ r}), where C̆ is
the reduced covariance measure.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `A`: Type of radii array
- `T`: Type of C function values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Fields
- `radii::A`: Distances at which the C function is evaluated
- `value::T`: C function values
- `processinformation::IP`: Information about the processes
- `estimationinformation::IE`: Information about the estimation procedure

# Mathematical Background
For a stationary process with power spectral density f(k), the C function is:
C(r) = ∫ f(k) W(r,k) dk
where W(r,k) is a spatial weighting function depending on the dimension.

# Examples
```julia
# Compute C function from spectral estimate
cf = c_function(spectrum, radii=0.1:0.1:2.0)

# Direct computation from data
cf = c_function(data, region, radii=0.1:0.1:2.0, nk=(32,32), kmax=(0.5,0.5), tapers=tapers)
```
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
getshortbaseestimatename(::Type{<:CFunction}) = "C"
getbaseestimatename(::Type{<:CFunction}) = "C function"
getargument(f::CFunction) = f.radii
getestimate(f::CFunction) = f.value

# Public API

"""
    c_function(data, region; kwargs...)

Compute spatial C function directly from data and region.
"""
function c_function(data, region; kwargs...)::CFunction
    return c_function(spatial_data(data, region); kwargs...)
end

"""
    c_function(data::SpatialData; radii, nk, kmax, wavenumber_radii, rotational_method, spectra_kwargs...)

Compute spatial C function from spatial data.

First computes the power spectral density, then transforms to the C function
via inverse Fourier transform with appropriate spatial weighting.

# Arguments
- `data::SpatialData`: Input spatial data
- `radii`: Distances at which to evaluate the C function
- `wavenumber_radii`: Radial wavenumbers for rotational averaging (default: from nk, kmax)
- `rotational_method`: Kernel for rotational averaging (default: from nk, kmax)
- `spectra_kwargs...`: Additional arguments passed to `spectra`

# Returns
A `CFunction` object containing the spatial C function.
"""
function c_function(data::SpatialData; radii = default_radii(data),
        rotational_wavenumber_radii = nothing,
        rotational_method = nothing, kwargs...)::CFunction
    spectrum = spectra(data; kwargs...)
    wavenumber_radii_processed = process_c_rotational_radii(
        spectrum, rotational_wavenumber_radii)
    rotational_method_processed = process_c_rotational_kernel(spectrum, rotational_method)
    return c_function(
        spectrum; radii = radii, wavenumber_radii = wavenumber_radii_processed,
        rotational_method = rotational_method_processed)
end
process_c_rotational_kernel(spectrum, ::Nothing) = default_c_rotational_kernel(spectrum)
process_c_rotational_kernel(spectrum, kernel) = kernel

process_c_rotational_radii(spectrum, ::Nothing) = default_c_rotational_radii(spectrum)
process_c_rotational_radii(spectrum, radii) = radii

default_c_rotational_kernel(spectrum) = NoRotational()
default_c_rotational_radii(spectrum) = default_rotational_radii(spectrum)

function default_radii(data::SpatialData)
    region = getregion(data)
    short_side = Meshes.ustrip(minimum(sides(region)))
    return range(0, short_side / 3, length = 100)
end

"""
    c_function(spectrum::Spectra; radii, wavenumber_radii, rotational_method)

Compute spatial C function from a spectral estimate.

# Arguments
- `spectrum::Spectra`: Input power spectral density estimate
- `radii`: Distances for C function evaluation
- `wavenumber_radii`: Radial wavenumbers for rotational averaging (default: from spectrum)
- `rotational_method`: Smoothing kernel for rotational averaging (default: from spectrum)

# Returns
A `CFunction` object with the C function values.
"""
function c_function(spectrum::Spectra; radii = default_radii(data),
        wavenumber_radii = process_c_rotational_kernel(spectrum, nothing),
        rotational_method = process_c_rotational_kernel(spectrum, nothing))
    rot_spec = rotational_estimate(
        spectrum, radii = wavenumber_radii, kernel = rotational_method)
    return _c_function(rot_spec, radii)
end

function c_function(spectrum::RotationalSpectra; radii = default_radii(data))
    _c_function(spectrum, radii)
end

function partial_c_function(data, region; kwargs...)::CFunction{PartialTrait}
    return partial_c_function(spatial_data(data, region); kwargs...)
end

function partial_c_function(
        data::SpatialData; radii = default_radii(data), spectra_kwargs...)::CFunction{PartialTrait}
    f_mt = partial_spectra(data; spectra_kwargs...)
    return c_function(f_mt; radii = radii)
end

function partial_c_function(
        spectrum::NormalOrRotationalSpectra{PartialTrait}; kwargs...
)::CFunction{PartialTrait}
    return c_function(spectrum; kwargs...)
end

function partial_c_function(
        spectrum::NormalOrRotationalSpectra{MarginalTrait}; kwargs...
)::CFunction{PartialTrait}
    return c_function(partial_spectra(spectrum); kwargs...)
end

# Internal implementation

function _c_function(spectrum::NormalOrRotationalSpectra{E}, radii) where {E}
    out = preallocate_c_output(spectrum, radii)
    store = precompute_c_weights(spectrum, radii)

    value = _sdf2C!(out, store, spectrum, radii)
    return CFunction{E}(
        radii, value, getprocessinformation(spectrum), getestimationinformation(spectrum))
end

# Core transformation functions

function _sdf2C!(out, store, spectrum::NormalOrRotationalSpectra, radii)
    power = getestimate(spectrum)
    zero_atom = getprocessinformation(spectrum).atoms
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

function preallocate_c_output(spectrum::NormalOrRotationalSpectra, radii)
    return _preallocate_c_output(getestimate(spectrum), process_trait(spectrum), radii)
end
function _preallocate_c_output(
        power, ::Union{SingleProcessTrait, MultipleTupleTrait}, radii)
    return zeros(real(eltype(power)), length(radii))
end
function _preallocate_c_output(power, ::MultipleVectorTrait, radii)
    return zeros(real(eltype(power)), size(power, 1), size(power, 2), length(radii))
end

function precompute_c_weights(spectra::Spectra{E, D}, radii::AbstractVector) where {E, D}
    wavenumber = getargument(spectra)
    power_size = size(getestimate(spectra))
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
    wavenumber = getargument(spectra)
    power_size = size(getestimate(spectra))
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

Compute spatial weighting function for D-dimensional c function.

These functions represent the Fourier transform of spherical/circular domains.
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

Compute weighting function for isotropic C functions.

For isotropic estimates, the weighting accounts for the radial integration
that has already been performed.
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
