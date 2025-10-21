"""
    c_function(spectrum::Spectra; radii, wavenumber_radii, rotational_method) -> CFunction
    c_function(spectrum::RotationalSpectra; radii) -> CFunction
    c_function(data, region; kwargs...) -> CFunction
    c_function(data::SpatialData; kwargs...) -> CFunction

Compute spatial C function from spectral estimates or directly from spatial data.

The C function is the reduced covariance measure evaluated on punctured balls as a
function of distance. For a stationary process at radius `r`, it represents
C(r) = C̆({x : 0 < ||x|| ≤ r}), where C̆ is the reduced covariance measure. The C function
is computed via inverse Fourier transform of the power spectral density with appropriate
spatial weighting functions that depend on the dimension.

# Arguments
- `spectrum::Spectra`: Input power spectral density estimate
- `spectrum::RotationalSpectra`: Rotationally averaged spectral estimate
- `data`: Spatial data for direct C function computation
- `region::Meshes.Geometry`: Spatial region for direct computation

# Keywords
- `radii::AbstractVector`: Distances at which to evaluate the C function (required)
- `wavenumber_radii::Union{Nothing, AbstractVector} = nothing`: Radial wavenumbers for
    rotational averaging. If `nothing`, defaults are computed from the spectrum grid.
- `rotational_method = nothing`: Smoothing kernel for rotational averaging. If `nothing`,
    uses `NoRotational()` (no smoothing).

When computing directly from data, all keywords from [`spectra`](@ref) are also supported:
- `nk`: Number of wavenumbers in each dimension
- `kmax`: Maximum wavenumber in each dimension
- `dk`: Wavenumber spacing in each dimension
- `tapers`: Taper functions to use
- `nw = 3`: Space-bandwidth product for taper generation
- `mean_method = DefaultMean()`: Method for mean estimation

# Returns
- `CFunction`: A C function estimate containing:
  - `radii`: Distance values where C function is evaluated
  - `value`: C function values (real-valued)
  - `processinformation`: Information about the analyzed processes
  - `estimationinformation`: Details about the estimation procedure

# Mathematical Details
For a stationary process with power spectral density f(k), the C function is computed as:

C(r) = ∫ f(k) W(r,k) dk

where W(r,k) is the spatial weighting function that depends on dimension D:
- **1D**: W(r,k) = 2r sinc(2r|k|)
- **2D**: W(r,k) = (r/|k|) J₁(2πr|k|) for |k| > 0, πr² for k = 0
- **3D**: W(r,k) = (r/|k|)^(3/2) J₃/₂(2πr|k|) for |k| > 0, (4π/3)r³ for k = 0

where J_ν is the Bessel function of the first kind of order ν, and sinc(x) = sin(πx)/(πx).

# Notes
- The C function measures spatial correlation as a function of distance
- Results are always real-valued due to the spatial integration
- For anisotropic spectra, rotational averaging is performed first
- Zero-atom corrections are automatically applied when available
- Use [`partial_c_function`](@ref) for partial C functions from partial spectra

# Examples
```julia
# Compute C function from existing spectral estimate
spec = spectra(data; kmax = 0.5, nw = 3)
radii = 0.1:0.1:2.0
cf = c_function(spec; radii = radii)

# Direct computation from data and region
cf = c_function(data, region; radii = 0.1:0.1:2.0, nk = (64, 64), kmax = (0.5, 0.5))

# With custom rotational averaging parameters
wavenumber_radii = 0.05:0.05:0.5
cf = c_function(spec; radii = radii, wavenumber_radii = wavenumber_radii)

# From rotational spectrum (no additional averaging needed)
rot_spec = rotational_spectra(data; kmax = 0.5, nw = 3)
cf = c_function(rot_spec; radii = radii)

# Access C function values
correlation_at_1km = cf.value[findfirst(r -> r ≈ 1.0, cf.radii)]
```

See also: [`spectra`](@ref), [`partial_c_function`](@ref), [`rotational_spectra`](@ref)
"""
c_function

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
getshortbaseestimatename(::Type{<:CFunction}) = "C"
getbaseestimatename(::Type{<:CFunction}) = "C function"
get_evaluation_points(f::CFunction) = f.radii
get_estimates(f::CFunction) = f.value

# Public API

function c_function(data, region; kwargs...)::CFunction
    return c_function(spatial_data(data, region); kwargs...)
end

function c_function(data::SpatialData; radii = nothing,
        rotational_wavenumber_radii = nothing,
        rotational_method = nothing, kwargs...)::CFunction
    spectrum = spectra(data; kwargs...)
    wavenumber_radii_processed = process_c_rotational_radii(
        spectrum, rotational_wavenumber_radii)
    rotational_method_processed = process_c_rotational_kernel(spectrum, rotational_method)
    radii = _validate_radii(radii, data)
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

function c_function(spectrum::Spectra; radii,
        wavenumber_radii = process_c_rotational_kernel(spectrum, nothing),
        rotational_method = process_c_rotational_kernel(spectrum, nothing))
    rot_spec = rotational_estimate(
        spectrum, radii = wavenumber_radii, kernel = rotational_method)
    return _c_function(rot_spec, radii)
end

function c_function(spectrum::RotationalSpectra; radii)
    _c_function(spectrum, radii)
end

"""
    partial_c_function(data, region; kwargs...) -> CFunction{PartialTrait}
    partial_c_function(data::SpatialData; kwargs...) -> CFunction{PartialTrait}
    partial_c_function(spectrum::NormalOrRotationalSpectra{PartialTrait}; kwargs...) -> CFunction{PartialTrait}
    partial_c_function(spectrum::NormalOrRotationalSpectra{MarginalTrait}; kwargs...) -> CFunction{PartialTrait}

Compute spatial partial C function from spectral estimates or directly from spatial data.

The partial C function represents the direct spatial correlation structure after removing
the linear influence of all other processes. It is computed from partial spectral estimates
using the same inverse Fourier transform approach as the regular C function, but applied
to partial spectra rather than marginal spectra.

# Arguments
- `data`: Spatial data for direct partial C function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `spectrum`: Spectral estimate (partial or marginal trait)

# Keywords
- `radii::AbstractVector`: Distances at which to evaluate the partial C function (required)
- Additional keywords: All keywords from [`partial_spectra`](@ref) and [`c_function`](@ref)

# Returns
- `CFunction{PartialTrait}`: A partial C function estimate with the same structure as
    regular C functions but representing direct spatial relationships.

# Notes
- Automatically converts marginal spectra to partial spectra when needed
- Uses the same spatial weighting functions as regular C functions
- Results represent direct spatial correlations after removing indirect effects
- For partial trait inputs, applies C function transformation directly

# Examples
```julia
# Direct computation from data
radii = 0.1:0.1:2.0
partial_cf = partial_c_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing partial spectral estimate
partial_spec = partial_spectra(data; kmax = 0.5, nw = 3)
partial_cf = partial_c_function(partial_spec; radii = radii)

# From marginal spectrum (automatically converts to partial)
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
partial_cf = partial_c_function(marginal_spec; radii = radii)
```

See also: [`c_function`](@ref), [`partial_spectra`](@ref)
"""
function partial_c_function(data, region; kwargs...)::CFunction{PartialTrait}
    return partial_c_function(spatial_data(data, region); kwargs...)
end

function partial_c_function(
        data::SpatialData; radii = nothing, spectra_kwargs...)::CFunction{PartialTrait}
    radii = _validate_radii(radii, data)
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
        radii, value, get_process_information(spectrum), get_estimation_information(spectrum))
end

# Core transformation functions

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

function preallocate_c_output(spectrum::NormalOrRotationalSpectra, radii)
    return _preallocate_c_output(get_estimates(spectrum), process_trait(spectrum), radii)
end
function _preallocate_c_output(
        power, ::Union{SingleProcessTrait, MultipleTupleTrait}, radii)
    return zeros(real(eltype(power)), length(radii))
end
function _preallocate_c_output(power, ::MultipleVectorTrait, radii)
    return zeros(real(eltype(power)), size(power, 1), size(power, 2), length(radii))
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
