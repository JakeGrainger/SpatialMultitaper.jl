struct RotationalEstimate{T, S, D, P} <: IsotropicEstimate{D, P}
    radii::T
    estimate::S
    function RotationalEstimate(a::AnisotropicEstimate{D, P}, radii, kernel) where {D, P}
        estimate = smoothed_rotational(getargument(a), getestimate(a), radii, kernel)
        new{typeof(radii), typeof(estimate), D, P}(radii, estimate)
    end
    function RotationalEstimate(
            radii::T, estimate::S, ::Val{D}, ::Val{P}) where {T, S, D, P}
        new{T, S, D, P}(radii, estimate)
    end
end
getargument(f::IsotropicEstimate) = f.radii
getestimate(f::IsotropicEstimate) = f.estimate
getextrafields(::RotationalEstimate{R, S, D, P}) where {R, S, D, P} = (Val{D}(), Val{P}())

struct NoRotational end # just indicates not to do rotational averaging
struct GaussKernel{T <: Real}
    bw::T
end
function (f::GaussKernel)(x)
    exp(-(x^2) / (2 * f.bw^2)) / f.bw
end
struct RectKernel{T <: Real}
    bw::T
end
function (f::RectKernel)(x)
    (abs(x / f.bw) < 1 / 2) / f.bw
end

_rotational_estimate(a, radii, ::NoRotational) = a
_rotational_estimate(a, radii, kernel) = RotationalEstimate(a, radii, kernel)
function rotational_estimate(
        f::AnisotropicEstimate{D, P};
        radii,
        kernel = Rect(mean(diff(radii)))
) where {D, P}
    _rotational_estimate(f, radii, kernel)
end

function smoothed_rotational(x, y, radii, kernel)
    @assert size(y) == length.(x)
    xitr = Iterators.ProductIterator(x)
    return [sum(f * kernel(norm(u) - r) for (u, f) in zip(xitr, y)) /
            sum(kernel(norm(u) - r) for u in xitr) for r in radii]
end

function default_rotational_radii(nfreq, fmax)
    default_rotational_radii(choose_freq_1d.(nfreq, fmax))
end

function default_rotational_radii(s::SpectralEstimate)
    return default_rotational_radii(getargument(s))
end

"""
    default_rotational_radii(freq::NTuple{D,AbstractVector{<:Real}}) where {D}
    default_rotational_radii(s::SpectralEstimate)
    default_rotational_radii(nfreq, fmax)

Constructs a default set of rotational averaging radii based on the frequency vectors.
The maximum radius is set to the minimum of the maximum frequencies in each dimension,
and the number of radii is set to the maximum length of the frequency vectors.
"""
function default_rotational_radii(freq::NTuple{D, AbstractVector{<:Real}}) where {D}
    max_freq = minimum(x -> maximum(abs, x), freq)
    n_freq = maximum(length, freq)
    zero_range = range(zero(eltype(max_freq)), stop = max_freq, length = n_freq)
    used_range = range(step(zero_range) / 2, stop = max_freq, step = step(zero_range)) # want to integrate with endpoints at zero range steps
    return used_range
end

"""
    default_rotational_kernel(freq_radii::OrdinalRange)

Constructs default rotational kernel based on the frequency radii.
The bandwidth is set to twice the step size of the frequency radii, and kernel is rectangular.
"""
function default_rotational_kernel(freq_radii::Union{StepRangeLen, OrdinalRange})
    bw = 2 * step(freq_radii) # ensures that kernels overlap with at least one evaluation point
    return RectKernel(bw)
end
