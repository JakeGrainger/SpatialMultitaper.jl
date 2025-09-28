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
cf = c_function(data, region, radii=0.1:0.1:2.0, nfreq=(32,32), fmax=(0.5,0.5), tapers=tapers)
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
    c_function(data::SpatialData; radii, nfreq, fmax, freq_radii, rotational_method, spectra_kwargs...)

Compute spatial C function from spatial data.

First computes the power spectral density, then transforms to the C function
via inverse Fourier transform with appropriate spatial weighting.

# Arguments
- `data::SpatialData`: Input spatial data
- `radii`: Distances at which to evaluate the C function
- `nfreq`: Number of frequencies per dimension for spectral estimation
- `fmax`: Maximum frequency per dimension for spectral estimation
- `freq_radii`: Radial frequencies for rotational averaging (default: from nfreq, fmax)
- `rotational_method`: Kernel for rotational averaging (default: from nfreq, fmax)
- `spectra_kwargs...`: Additional arguments passed to spectral estimation

# Returns
A `CFunction` object containing the spatial C function.
"""
function c_function(
        data::SpatialData; radii, nfreq, fmax,
        freq_radii = default_rotational_radii(nfreq, fmax),
        rotational_method = default_rotational_kernel(nfreq, fmax),
        spectra_kwargs...
)::CFunction
    spectrum = spectra(data; nfreq, fmax, spectra_kwargs...)
    return c_function(spectrum; radii = radii, freq_radii = freq_radii,
        rotational_method = rotational_method)
end

"""
    c_function(spectrum::Spectra; radii, freq_radii, rotational_method)

Compute spatial C function from a spectral estimate.

# Arguments
- `spectrum::Spectra`: Input power spectral density estimate
- `radii`: Distances for C function evaluation
- `freq_radii`: Radial frequencies for rotational averaging (default: from spectrum)
- `rotational_method`: Smoothing kernel for rotational averaging (default: from spectrum)

# Returns
A `CFunction` object with the C function values.
"""
function c_function(
        spectrum::Spectra; radii,
        freq_radii = default_rotational_radii(spectrum),
        rotational_method = default_rotational_kernel(spectrum)
)::CFunction
    return _c_function(spectrum, radii, freq_radii, rotational_method)
end

"""
    c_function(spectrum::RotationalSpectra{E}; radii) where {E}

Compute C function from an already rotationally averaged spectrum.
"""
function c_function(spectrum::RotationalSpectra{E}; radii)::CFunction{E} where {E}
    value = sdf2C(spectrum, radii)
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
    partial_c_function(data::SpatialData; radii, nfreq, fmax, spectra_kwargs...)

Compute partial spatial C function from spatial data.

Partial C functions show direct spatial relationships with the influence
of other processes removed.
"""
function partial_c_function(
        data::SpatialData; radii, nfreq, fmax, spectra_kwargs...
)::CFunction{PartialTrait}
    f_mt = partial_spectra(data; nfreq, fmax, spectra_kwargs...)
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
    value = sdf2C(rot_spec, radii)
    return CFunction{E}(
        radii, value, getprocessinformation(spectrum), getestimationinformation(spectrum))
end

# Core transformation functions

"""
    sdf2C(f, radii::AbstractVector{<:Number})

Convert spectral density to C function at multiple radii.

This function computes an inverse transform of the power spectral density
to obtain the spatial C function.

# Arguments
- `f`: Spectral estimate (anisotropic or isotropic)
- `radii`: Vector of distances at which to evaluate the C function

# Returns
Vector of C function values corresponding to the input radii.
"""
function sdf2C(f, radii::AbstractVector{<:Number})
    return [_sdf2C(f, radius) for radius in radii]
end

"""
    _sdf2C(f::Spectra, radius::Number)

Convert spectral density to C function at a single radius.

Uses appropriate weighting functions based on spatial dimension and handles
zero-atom corrections when present.
"""
function _sdf2C(f::Spectra, radius::Number)
    freq = getargument(f)
    spectra = getestimate(f)
    zero_atom = getprocessinformation(f).atoms
    return _sdf2C_anisotropic(freq, spectra, process_trait(f), zero_atom, radius)
end

"""
    _sdf2C(f::IsotropicEstimate{E, D, P}, radius::Number) where {E, D, P}

Convert isotropic spectral density to C function at a single radius.

For isotropic estimates, uses specialized weighting functions that account for
the radial symmetry.
"""
function _sdf2C(f::IsotropicEstimate{E, D, P}, radius::Number) where {E, D, P}
    freq = getargument(f)
    spectra = getestimate(f)
    zero_atom = getprocessinformation(f).atoms
    return _sdf2C_isotropic(freq, spectra, process_trait(f), zero_atom, radius, Val{D}())
end

# Anisotropic correlation function computation

function _sdf2C_anisotropic(freq, power, ::Union{SingleProcessTrait, MultipleTupleTrait},
        zero_atom, radius::Number)
    frequency_spacing = prod(step, freq)

    c = sum(_compute_c_term(s, zero_atom, radius, k)
    for (s, k) in zip(power, Iterators.product(freq...)))

    return frequency_spacing * real(c)
end

function _sdf2C_anisotropic(freq, power, ::MultipleVectorTrait, zero_atom, radius::Number)
    if length(freq) > ndims(power)
        throw(DimensionMismatch("Frequency dimensions ($(length(freq))) cannot exceed power array dimensions ($(ndims(power)))"))
    end

    out = mapslices(
        z -> _sdf2C_anisotropic(freq, z, SingleProcessTrait(), zero_atom, radius),
        power; dims = (N - D + 1):ndims(power))
    return reshape(out, size(out)[1:(N - D)]) # The D spatial dimensions have been collapsed into one singleton (at one radius)
end

## _sdf2C_isotropic
function _sdf2C_isotropic(freq, power, ::Union{SingleProcessTrait, MultipleTupleTrait},
        zero_atom, radius::Number, ::Val{D}) where {D}
    spacing = step(freq)

    c = sum(_compute_c_term(s, zero_atom, radius, k, spacing, Val{D}())
    for (s, k) in zip(power, freq))

    return real(c)
end

function _sdf2C_isotropic(
        freq, power, ::MultipleVectorTrait, zero_atom, radius::Number, ::Val{D}) where {D}
    if length(freq) > ndims(power)
        throw(DimensionMismatch("Frequency dimensions ($(length(freq))) cannot exceed power array dimensions ($(ndims(power)))"))
    end

    out = mapslices(
        z -> _sdf2C_isotropic(freq, z, SingleProcessTrait(), zero_atom, radius, Val{D}()),
        power; dims = (N - D + 1):ndims(power))
    return reshape(out, size(out)[1:(N - D)]) # The D spatial dimensions have been collapsed into one singleton (at one radius)
end

# Helper functions for C function computation

"""
    _compute_c_term(s, zero_atom, radius, k, spacing)

Compute individual terms in the C function sum.
"""
function _compute_c_term(s, zero_atom, radius, k)
    weight = _sphere_weight(radius, k, Val{length(k)}())
    return (s - zero_atom) * weight
end
function _compute_c_term(s, ::Nothing, radius, k)
    weight = _sphere_weight(radius, k, Val{length(k)}())
    return s * weight
end

function _compute_c_term(s, zero_atom, radius, k, spacing, ::Val{D}) where {D}
    weight = _isotropic_weight(radius, k, spacing, Val{D}())
    return (s - zero_atom) * weight
end
function _compute_c_term(s, ::Nothing, radius, k, spacing, ::Val{D}) where {D}
    weight = _isotropic_weight(radius, k, spacing, Val{D}())
    return s * weight
end

# Spatial weighting functions

"""
    _sphere_weight(r, u, ::Val{D})

Compute spatial weighting function for D-dimensional c function.

These functions represent the Fourier transform of spherical/circular domains.
"""
function _sphere_weight(r, u, ::Val{1})
    x = norm(u)
    return 2r * sinc(2r * x)
end

function _sphere_weight(r, u, ::Val{2})
    x = norm(u)
    return (x < 1e-10) ? (π * r^2) : ((r / x) * besselj1(2π * r * x))
end

function _sphere_weight(r, u, ::Val{D}) where {D}
    x = norm(u)
    if x < 1e-10
        # Handle singularity at origin using ball measure
        return unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), r))
    else
        return (r / x)^(D / 2) * besselj(D / 2, 2π * r * x)
    end
end

"""
    _isotropic_weight(r, k, spacing, ::Val{D})

Compute weighting function for isotropic C functions.

For isotropic estimates, the weighting accounts for the radial integration
that has already been performed.
"""
function _isotropic_weight(r, k, spacing, ::Val{2})
    half_spacing = spacing / 2
    return besselj0(2π * r * (k - half_spacing)) - besselj0(2π * r * (k + half_spacing))
end

function _isotropic_weight(r, k, spacing, ::Val{D}) where {D}
    error("Isotropic weighting not implemented for D != 2")
end
