"""
    RotationalEstimate{E, D, S, A, T, IP, IE} <: IsotropicEstimate{E, D}

An estimate that has undergone rotational averaging to produce isotropic results.

Rotational averaging transforms anisotropic estimates (which vary with direction)
into isotropic estimates that depend only on radial distance from the origin.
This is commonly used in spectral analysis to study scale-dependent properties
independent of orientation.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `S`: Type of the original anisotropic estimate
- `A`: Type of radii array
- `T`: Type of the rotationally averaged estimate values
- `IP`: Type of process information
- `IE`: Type of estimation information

# Fields
- `radii::A`: Radial distances at which the estimate is evaluated
- `estimate::T`: Rotationally averaged estimate values
- `processinformation::IP`: Information about the processes
- `estimationinformation::IE`: Information about the estimation procedure

# Examples
```julia
# Create rotational estimate with default parameters
rot_spec = rotational_estimate(anisotropic_spectrum)

# Create with custom radii and kernel
radii = 0.1:0.1:1.0
kernel = GaussKernel(0.05)
rot_spec = rotational_estimate(spectrum, radii=radii, kernel=kernel)
```
"""
struct RotationalEstimate{E, D, S, A, T, IP, IE} <: IsotropicEstimate{E, D}
    radii::A
    estimate::T
    processinformation::IP
    estimationinformation::IE
    function RotationalEstimate{E, S}(
            radii::A, estimate::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, S, A, T, D, IE}
        checkinputs(radii, estimate, processinfo)
        IP = typeof(processinfo)
        return new{E, D, S, A, T, IP, IE}(
            radii, estimate, processinfo, estimationinfo)
    end
end

"""
    getargument(f::RotationalEstimate)

Get the radii at which the rotational estimate is evaluated.
"""
getargument(f::RotationalEstimate) = f.radii

"""
    getestimate(f::RotationalEstimate)

Get the rotationally averaged estimate values.
"""
getestimate(f::RotationalEstimate) = f.estimate

"""
    getoriginaltype(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S}

Extract the original estimate type before rotational averaging.
"""
getoriginaltype(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S} = S

"""
    rotational_estimate(est::AnisotropicEstimate{E}; radii, kernel) where {E}

Compute a rotationally averaged estimate from an anisotropic estimate.

Rotational averaging integrates the estimate over circular annuli at different
radii, producing an isotropic result that depends only on distance from the origin.
This is useful for analyzing scale-dependent properties without directional bias.

# Arguments
- `est::AnisotropicEstimate`: The anisotropic estimate to average
- `radii`: Radial distances for evaluation (default: `default_rotational_radii(est)`)
- `kernel`: Smoothing kernel for averaging (default: `default_rotational_kernel(est)`)

# Returns
A `RotationalEstimate` containing the rotationally averaged values.
The exception is if kernel is `NoRotational`, in which case the input estimate is returned.
This is used to skip rotational averaging in some downstream cases.

# Examples
```julia
# Basic rotational averaging
rot_est = rotational_estimate(spectrum)

# With custom parameters
radii = range(0.1, 2.0, length=50)
kernel = GaussKernel(0.1)  # Gaussian smoothing with bandwidth 0.1
rot_est = rotational_estimate(spectrum, radii=radii, kernel=kernel)
```
"""
function rotational_estimate(
        est::AnisotropicEstimate; radii = default_rotational_radii(est),
        kernel = default_rotational_kernel(est))
    return _rotational_estimate(est, radii, kernel)
end

"""
    getestimatename(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S}

Get the name of the estimate, reflecting its rotational and original traits.

The name indicates the order of application of rotational and partial traits.
"""
function getestimatename(::Type{<:RotationalEstimate{E, D, S}}) where {E, D, S}
    original_name = getbaseestimatename(S)

    if E === MarginalTrait
        return "rotational " * original_name
    elseif E === PartialTrait
        # Check if original was marginal or partial to determine order
        original_trait = _extract_estimate_trait(S)
        if original_trait === MarginalTrait
            return "partial rotational " * original_name  # Partial applied after rotational
        else
            return "rotational partial " * original_name  # Rotational applied after partial
        end
    end
end
# Helper function to extract trait from type
_extract_estimate_trait(::Type{<:AbstractEstimate{T}}) where {T} = T

"""
    _construct_estimate_subset(
        ::Type{<:RotationalEstimate{E, D, S}},
        trait::Type{<:EstimateTrait},
        argument, estimate, processinfo, estimationinfo
    ) where {E, D, S}

Construct a subset estimate with a different trait from an existing RotationalEstimate.

This is used to create partial or marginal estimates from a rotational estimate.
"""
function _construct_estimate_subset(
        ::Type{<:RotationalEstimate{E, D, S}},
        trait::Type{<:EstimateTrait},
        argument, estimate, processinfo, estimationinfo
) where {E, D, S}
    # For RotationalEstimate, we need {E, S} constructor
    return RotationalEstimate{trait, S}(
        argument, estimate, processinfo, estimationinfo
    )
end

## Smoothing kernels for rotational averaging

"""
    NoRotational

Indicates that no rotational averaging should be performed.
"""
struct NoRotational end # just indicates not to do rotational averaging

"""
    GaussKernel{T <: Real}

Gaussian smoothing kernel for rotational averaging.

# Fields
- `bw::T`: Bandwidth parameter controlling the width of the Gaussian

# Mathematical Form
K(x) = exp(-x²/(2σ²)) / σ where σ is the bandwidth.
"""
struct GaussKernel{T <: Real}
    bw::T
end

"""
    (f::GaussKernel)(x)

Evaluate the Gaussian kernel at distance x.
"""
function (f::GaussKernel)(x)
    return exp(-(x^2) / (2 * f.bw^2)) / f.bw
end

"""
    RectKernel{T <: Real}

Rectangular (uniform) smoothing kernel for rotational averaging.

# Fields
- `bw::T`: Bandwidth parameter controlling the width of the rectangular window

# Mathematical Form
K(x) = 1/h if |x/h| < 1/2, 0 otherwise, where h is the bandwidth.
"""
struct RectKernel{T <: Real}
    bw::T
end

"""
    (f::RectKernel)(x)

Evaluate the rectangular kernel at distance x.
"""
function (f::RectKernel)(x)
    return (abs(x / f.bw) < 1 / 2) / f.bw
end

# Internal implementation functions

"""
    _rotational_estimate(a::AbstractEstimate, radii, ::NoRotational)

Pass-through function when no rotational averaging is requested.
"""
_rotational_estimate(a::AbstractEstimate, radii, ::NoRotational) = a

"""
    _rotational_estimate(a::AbstractEstimate, radii, kernel)

Internal function to compute rotational averaging with a given kernel.
"""
function _rotational_estimate(est::AbstractEstimate{E}, radii, kernel) where {E}
    rot_est = _smoothed_rotational(
        getargument(est), getestimate(est), process_trait(est), radii, kernel)
    processinfo = getprocessinformation(est)
    estimationinfo = getestimationinformation(est)
    return RotationalEstimate{E, typeof(est)}(radii, rot_est, processinfo, estimationinfo)
end

"""
    _smoothed_rotational(x::NTuple{D}, y::AbstractArray{T, D}, radii, kernel) where {D, T}

Core rotational averaging computation.

Computes weighted averages over circular annuli using the specified kernel.
"""
function _smoothed_rotational(x::NTuple{D}, y::AbstractArray{T, D},
        ::Union{SingleProcessTrait, MultipleTupleTrait}, radii,
        kernel) where {D, T <: Union{<:Number, <:SMatrix}}
    @argcheck length(x) == ndims(y)
    @argcheck size(y) == length.(x)
    xitr = Iterators.ProductIterator(x)
    return [sum(f * kernel(norm(u) - r) for (u, f) in zip(xitr, y)) /
            sum(kernel(norm(u) - r) for u in xitr) for r in radii]
end

function _smoothed_rotational(
        x::NTuple{D}, y::AbstractArray{<:Number, N}, ::MultipleVectorTrait, radii, kernel) where {
        D, N}
    @argcheck length(x) <= ndims(y)
    # assumes the that last D dimensions are the spatial ones
    @argcheck size(y)[3:end] == length.(x)
    out = mapslices(
        z -> _smoothed_rotational(x, z, SingleProcessTrait(), radii, kernel), y; dims = (N + 1 - D):ndims(y))
    return reshape(out, size(out)[1:(N + 1 - D)]) # The D spatial dimensions have been collapsed into one
end

"""
    default_rotational_radii(wavenumber::NTuple{D,AbstractVector{<:Real}}) where {D}
    default_rotational_radii(s::AbstractEstimate)
    default_rotational_radii(nk, kmax)

Construct default radii for rotational averaging based on wavenumber vectors.

The default radii are chosen to span from near zero to the minimum of the
maximum wavenumbers across dimensions, with the number of points equal to
the maximum wavenumber vector length.

# Arguments
- `wavenumber`: Tuple of wavenumber vectors for each dimension
- `s`: An estimate from which to extract wavenumber information
- `nk`, `kmax`: Number of wavenumbers and maximum wavenumber per dimension

# Returns
A range of radii suitable for rotational averaging.

# Algorithm
- Maximum radius = min(max_wavenumber_per_dimension)
- Number of radii = max(length_per_dimension)
- Radii are offset by half a step to avoid the origin singularity
"""
function default_rotational_radii(nk, kmax)
    return default_rotational_radii(_choose_wavenumbers_1d.(nk, kmax))
end

function default_rotational_radii(s::AbstractEstimate)
    return default_rotational_radii(getargument(s))
end

function default_rotational_radii(wavenumber::NTuple{D, AbstractVector{<:Real}}) where {D}
    max_wavenumber = minimum(x -> maximum(abs, x), wavenumber)
    n_wavenumber = maximum(length, wavenumber)
    zero_range = range(
        zero(eltype(max_wavenumber)), stop = max_wavenumber, length = n_wavenumber)
    # Offset by half step to integrate with endpoints at zero range steps
    used_range = range(step(zero_range) / 2, stop = max_wavenumber, step = step(zero_range))
    return used_range
end

"""
    default_rotational_kernel(est::AbstractEstimate)
    default_rotational_kernel(nk, kmax)
    default_rotational_kernel(wavenumber::NTuple{D, AbstractVector{<:Real}}) where {D}

Construct a default smoothing kernel for rotational averaging.

The default kernel is a rectangular kernel with bandwidth twice the maximum
wavenumber step size across dimensions. This provides reasonable smoothing
while preserving wavenumber resolution.

# Returns
A `RectKernel` with bandwidth = 2 * max(wavenumber_steps)
"""
function default_rotational_kernel(est::AbstractEstimate)
    return default_rotational_kernel(getargument(est))
end

function default_rotational_kernel(nk, kmax)
    return default_rotational_kernel(_choose_wavenumbers_1d.(nk, kmax))
end

function default_rotational_kernel(wavenumber::NTuple{D, AbstractVector{<:Real}}) where {D}
    max_step = maximum(step, wavenumber)
    return RectKernel(2 * max_step)
end

function default_rotational_kernel(wavenumber::AbstractVector{<:Real})
    max_step = step(wavenumber)
    return RectKernel(2 * max_step)
end
