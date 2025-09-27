struct RotationalEstimate{E, D, P, Q, S, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}
    radii::A
    estimate::T
    processinformation::IP
    estimationinformation::IE
    function RotationalEstimate{E, S}(
            radii::A, estimate::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, S, A, T, D, IE}
        P, Q = checkinputs(radii, estimate, processinfo)
        IP = typeof(processinfo)
        new{E, D, P, Q, S, A, T, IP, IE}(radii, estimate, processinfo, estimationinfo)
    end
end

getargument(f::RotationalEstimate) = f.radii
getestimate(f::RotationalEstimate) = f.estimate

# Helper function for dispatch on original type
getoriginaltype(::Type{<:RotationalEstimate{E, D, P, Q, S}}) where {E, D, P, Q, S} = S

function rotational_estimate(
        est::AnisotropicEstimate{E}; radii = default_rotational_radii(est),
        kernel = default_rotational_kernel(est)) where {E}
    rot_est = _rotational_estimate(est, radii, kernel)
    processinfo = getprocessinformation(est)
    estimationinfo = getestimationinformation(est)
    return RotationalEstimate{E, typeof(est)}(radii, rot_est, processinfo, estimationinfo)
end

function getestimatename(::Type{<:RotationalEstimate{E, D, P, Q, S}}) where {E, D, P, Q, S}
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

function _construct_estimate_subset(
        ::Type{<:RotationalEstimate{E, D, P, Q, S}},
        trait::Type{<:EstimateTrait},
        argument, estimate, processinfo, estimationinfo
) where {E, D, P, Q, S}
    # For RotationalEstimate, we need {E, S} constructor
    return RotationalEstimate{trait, S}(
        argument, estimate, processinfo, estimationinfo
    )
end

## internals

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

_rotational_estimate(a::AbstractEstimate, radii, ::NoRotational) = a
function _rotational_estimate(a::AbstractEstimate, radii, kernel)
    return smoothed_rotational(getargument(a), getestimate(a), radii, kernel)
end

function _smoothed_rotational(
        x::NTuple{D}, y::AbstractArray{T, D}, radii, kernel) where {
        D, T <: Union{<:Number, <:SMatrix}}
    xitr = Iterators.ProductIterator(x)
    return [sum(f * kernel(norm(u) - r) for (u, f) in zip(xitr, y)) /
            sum(kernel(norm(u) - r) for u in xitr) for r in radii]
end

function smoothed_rotational(
        x::NTuple{D}, y::AbstractArray{<:Number, D}, radii, kernel) where {D}
    @argcheck length(x) == ndims(y)
    @argcheck size(y) == length.(x)
    _smoothed_rotational(x, y, radii, kernel)
end
function smoothed_rotational(
        x::NTuple{D}, y::AbstractArray{<:SMatrix, D}, radii, kernel) where {D}
    @argcheck length(x) == ndims(y)
    @argcheck size(y) == length.(x)
    _smoothed_rotational(x, y, radii, kernel)
end
function smoothed_rotational(
        x::NTuple{D}, y::AbstractArray{<:Number, N}, radii, kernel) where {D, N}
    @argcheck length(x) == ndims(y) - 2
    @argcheck size(y)[3:end] == length.(x)
    out = mapslices(z -> _smoothed_rotational(x, z, radii, kernel), y; dims = 3:ndims(y))
    return reshape(out, size(out)[1:(ndims(out) - 1)])
end

"""
    default_rotational_radii(freq::NTuple{D,AbstractVector{<:Real}}) where {D}
    default_rotational_radii(s::AbstractEstimate)
    default_rotational_radii(nfreq, fmax)

Constructs a default set of rotational averaging radii based on the frequency vectors.
The maximum radius is set to the minimum of the maximum frequencies in each dimension,
and the number of radii is set to the maximum length of the frequency vectors.
"""
function default_rotational_radii(nfreq, fmax)
    default_rotational_radii(_choose_frequencies_1d.(nfreq, fmax))
end

function default_rotational_radii(s::AbstractEstimate)
    return default_rotational_radii(getargument(s))
end
function default_rotational_radii(freq::NTuple{D, AbstractVector{<:Real}}) where {D}
    max_freq = minimum(x -> maximum(abs, x), freq)
    n_freq = maximum(length, freq)
    zero_range = range(zero(eltype(max_freq)), stop = max_freq, length = n_freq)
    used_range = range(step(zero_range) / 2, stop = max_freq, step = step(zero_range)) # want to integrate with endpoints at zero range steps
    return used_range
end

function default_rotational_kernel(est::AbstractEstimate)
    return default_rotational_kernel(getargument(est))
end
function default_rotational_kernel(nfreq, fmax)
    return default_rotational_kernel(_choose_frequencies_1d.(nfreq, fmax))
end
function default_rotational_kernel(freq::NTuple{D, AbstractVector{<:Real}}) where {D}
    max_step = maximum(step, freq)
    return RectKernel(2 * max_step)
end
