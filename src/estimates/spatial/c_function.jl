"""
    CFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}

Spatial C function estimate derived from spectral estimates.

The C function is the reduced covariance measure of a process evaluated on a punctured ball,
of a function of distance. So at radius `r` it is C(r) = C̆({x : 0 < ||x|| ≤ r}), where C̆ is
the reduced covariance measure.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `P`, `Q`: Process dimensions
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
struct CFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function CFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        P, Q = checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        return new{E, D, P, Q, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end

getbaseestimatename(::Type{<:CFunction}) = "C function"

"""
    getargument(f::CFunction)

Get the radii at which the C function is evaluated.
"""
getargument(f::CFunction) = f.radii

"""
    getestimate(f::CFunction)

Get the C function values.
"""
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
    c_function(data::SpatialData; radii, nk, kmax, freq_radii, rotational_method, spectra_kwargs...)

Compute spatial C function from spatial data.

First computes the power spectral density, then transforms to the C function
via inverse Fourier transform with appropriate spatial weighting.

# Arguments
- `data::SpatialData`: Input spatial data
- `radii`: Distances at which to evaluate the C function
- `nk`: Number of wavenumbers per dimension for spectral estimation
- `kmax`: Maximum wavenumber per dimension for spectral estimation
- `freq_radii`: Radial wavenumbers for rotational averaging (default: from nk, kmax)
- `rotational_method`: Kernel for rotational averaging (default: from nk, kmax)
- `spectra_kwargs...`: Additional arguments passed to spectral estimation

# Returns
A `CFunction` object containing the spatial C function.
"""
function c_function(data::SpatialData; radii, rotational_freq_radii = nothing,
        rotational_method = nothing, kwargs...)::CFunction
    spectrum = spectra(data; kwargs...)
    freq_radii_processed = process_c_rotational_radii(rotational_freq_radii; kwargs...)
    rotational_method_processed = process_c_rotational_kernel(rotational_method; kwargs...)
    return c_function(spectrum; radii = radii, freq_radii = freq_radii_processed,
        rotational_method = rotational_method_processed)
end
process_c_rotational_kernel(::Nothing; kwargs...) = default_c_rotational_kernel(; kwargs...)
process_c_rotational_kernel(kernel; kwargs...) = kernel
process_c_rotational_radii(::Nothing; kwargs...) = default_c_rotational_radii(; kwargs...)
process_c_rotational_radii(radii; kwargs...) = radii

default_c_rotational_kernel(args...; kwargs...) = NoRotational()
default_c_rotational_radii(spectrum) = default_rotational_radii(spectrum)
default_c_rotational_radii(; nk, kmax, kwargs...) = default_rotational_radii(nk, kmax)

"""
    c_function(spectrum::Spectra; radii, freq_radii, rotational_method)

Compute spatial C function from a spectral estimate.

# Arguments
- `spectrum::Spectra`: Input power spectral density estimate
- `radii`: Distances for C function evaluation
- `freq_radii`: Radial wavenumbers for rotational averaging (default: from spectrum)
- `rotational_method`: Smoothing kernel for rotational averaging (default: from spectrum)

# Returns
A `CFunction` object with the C function values.
"""
function c_function(
        spectrum::Spectra; radii,
        freq_radii = default_c_rotational_radii(spectrum),
        rotational_method = default_c_rotational_kernel(spectrum)
)::CFunction
    return _c_function(spectrum, radii, freq_radii, rotational_method)
end

"""
    c_function(spectrum::RotationalSpectra{E}; radii) where {E}

Compute C function from an already rotationally averaged spectrum.
"""
function c_function(spectrum::RotationalSpectra{E}; radii)::CFunction{E} where {E}
    value = _sdf2C(spectrum, radii)
    return CFunction{E}(
        radii, value, getprocessinformation(spectrum), getestimationinformation(spectrum))
end

# Partial correlation functions

"""
    partial_c_function(data, region; kwargs...)

Compute partial spatial C function directly from data and region.
"""
function partial_c_function(data, region; kwargs...)::CFunction{PartialTrait}
    return partial_c_function(spatial_data(data, region); kwargs...)
end

"""
    partial_c_function(data::SpatialData; radii, nk, kmax, spectra_kwargs...)

Compute partial spatial C function from spatial data.

Partial C functions show direct spatial relationships with the influence
of other processes removed.
"""
function partial_c_function(
        data::SpatialData; radii, spectra_kwargs...)::CFunction{PartialTrait}
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

"""
    _c_function(spectrum::Spectra{E}, radii, freq_radii, rotational_method) where {E}

Internal function to compute C function with rotational averaging.
"""
function _c_function(spectrum::Spectra{E}, radii, freq_radii, rotational_method) where {E}
    # Perform rotational averaging if needed
    rot_spec = rotational_estimate(spectrum, radii = freq_radii, kernel = rotational_method)
    value = _sdf2C(rot_spec, radii)
    return CFunction{E}(
        radii, value, getprocessinformation(spectrum), getestimationinformation(spectrum))
end

# Core transformation functions

"""
    _sdf2C(f::Spectra, radii::Number)

Convert spectral density to C function at some collection of radii.

Uses appropriate weighting functions based on spatial dimension and handles
zero-atom corrections when present.
"""
function _sdf2C(f::Spectra, radii)
    wavenumber = getargument(f)
    spectra = getestimate(f)
    zero_atom = getprocessinformation(f).atoms
    return _sdf2C_anisotropic(wavenumber, spectra, process_trait(f), zero_atom, radii)
end

"""
    _sdf2C(f::IsotropicEstimate{E, D, P}, radius) where {E, D, P}

Convert isotropic spectral density to C function at some radii.

For isotropic estimates, weights are integrated over each interval in the integral being
approximated.
"""
function _sdf2C(f::IsotropicEstimate{E, D, P}, radii) where {E, D, P}
    wavenumber = getargument(f)
    spectra = getestimate(f)
    zero_atom = getprocessinformation(f).atoms
    return _sdf2C_isotropic(
        wavenumber, spectra, process_trait(f), zero_atom, radii, Val{D}())
end

# Anisotropic correlation function computation

function _sdf2C_anisotropic(
        wavenumber, power, ::Union{SingleProcessTrait, MultipleTupleTrait},
        zero_atom, radius::Number)
    wavenumber_spacing = prod(step, wavenumber)

    c = sum(_compute_c_term(s, zero_atom, radius, k)
    for (s, k) in zip(power, Iterators.product(wavenumber...)))

    return wavenumber_spacing * real(c)
end

function _sdf2C_anisotropic(
        wavenumber, power, trait::Union{SingleProcessTrait, MultipleTupleTrait},
        zero_atom, radii::AbstractVector)
    out = zeros(eltype(power), length(radii))
    for (i, radius) in enumerate(radii)
        out[i] = _sdf2C_anisotropic(wavenumber, power, trait, zero_atom, radius)
    end
    return out
end

function _sdf2C_anisotropic(
        wavenumber::NTuple{D}, power::AbstractArray{<:Number, N}, ::MultipleVectorTrait,
        zero_atom, radii::AbstractVector) where {D, N}
    if length(wavenumber) > ndims(power)
        throw(DimensionMismatch("Wavenumber dimensions ($(length(wavenumber))) cannot exceed power array dimensions ($(ndims(power)))"))
    end
    if !isnothing(zero_atom)
        @argcheck ndims(zero_atom) + length(wavenumber) == ndims(power)
    end

    out = zeros(eltype(power), size(power)[1:((N - D))]..., length(radii))
    for idx in CartesianIndices(size(power)[1:(N - D)])
        power_slice = view(power, idx, ntuple(Returns(:), D)...)
        zero_atom_slice = _slice_zero_atom(zero_atom, idx)
        for (i, radius) in enumerate(radii)
            out[idx, i] = _sdf2C_anisotropic(
                wavenumber, power_slice, SingleProcessTrait(), zero_atom_slice, radius)
        end
    end
    return out
end

## _sdf2C_isotropic
function _sdf2C_isotropic(
        wavenumber, power, trait::Union{SingleProcessTrait, MultipleTupleTrait},
        zero_atom, radii::AbstractVector, ::Val{D}) where {D}
    out = zeros(eltype(power), length(radii))
    for (i, radius) in enumerate(radii)
        out[i] = _sdf2C_isotropic(wavenumber, power, trait, zero_atom, radius, Val{D}())
    end
    return out
end

function _sdf2C_isotropic(
        wavenumber, power, ::Union{SingleProcessTrait, MultipleTupleTrait},
        zero_atom, radius::Number, ::Val{D}) where {D}
    spacing = step(wavenumber)

    c = sum(_compute_c_term(s, zero_atom, radius, k, spacing, Val{D}())
    for (s, k) in zip(power, wavenumber))

    return real(c)
end

function _sdf2C_isotropic(
        wavenumber, power::AbstractArray{<:Number, N}, ::MultipleVectorTrait,
        zero_atom, radii::AbstractVector, ::Val{D}) where {D, N}
    if !isnothing(zero_atom)
        @argcheck ndims(zero_atom) + 1 == ndims(power)
    end
    @argcheck size(power, ndims(power)) == length(wavenumber)

    out = zeros(eltype(power), size(power)[1:(N - 1)]..., length(radii))
    for idx in CartesianIndices(size(power)[1:(N - 1)])
        power_slice = view(power, idx, ntuple(Returns(:), D)...)
        zero_atom_slice = _slice_zero_atom(zero_atom, idx)
        for (i, radius) in enumerate(radii)
            out[idx, i] = _sdf2C_isotropic(
                wavenumber, power_slice, SingleProcessTrait(), zero_atom_slice, radius, Val{D}())
        end
    end
    return out
end

_slice_zero_atom(zero_atom, idx) = zero_atom[idx]
_slice_zero_atom(::Nothing, _) = nothing

# Helper functions for C function computation

"""
    _compute_c_term(s, zero_atom, radius, k, spacing)

Compute individual terms in the C function sum.
"""
function _compute_c_term(s, zero_atom, radius, k)
    weight = _anisotropic_c_weight(radius, k, Val{length(k)}())
    return (s - zero_atom) * weight
end
function _compute_c_term(s, ::Nothing, radius, k)
    weight = _anisotropic_c_weight(radius, k, Val{length(k)}())
    return s * weight
end

function _compute_c_term(s, zero_atom, radius, k, spacing, ::Val{D}) where {D}
    weight = _isotropic_c_weight(radius, k, spacing, Val{D}())
    return (s - zero_atom) * weight
end
function _compute_c_term(s, ::Nothing, radius, k, spacing, ::Val{D}) where {D}
    weight = _isotropic_c_weight(radius, k, spacing, Val{D}())
    return s * weight
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
