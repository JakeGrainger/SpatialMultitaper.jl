"""
    MarginallyTransformedEstimate{E, D, P, Q, N, S, A, T, IP, IE, F} <: AbstractEstimate{E, D, P, Q, N}

A wrapper for estimates that have undergone element-wise transformations.

This structure represents an estimate where a mathematical function has been applied
element-wise to the estimate values while preserving the argument structure and metadata.
Common examples include taking the real part, magnitude, phase, or logarithm of estimates.

# Type Parameters
- `E`: Estimate trait (e.g., `MarginalTrait`, `PartialTrait`)
- `D`: Spatial dimension
- `P`, `Q`: Process dimensions
- `N`: Number of wavenumber dimensions
- `S`: Type of the original estimate before transformation
- `A`: Type of the argument (e.g., wavenumber grid)
- `T`: Type of the transformed estimate values
- `IP`: Type of process information
- `IE`: Type of estimation information
- `F`: Type of the transformation function

# Fields
- `argument`: The argument (typically wavenumber) of the estimate
- `estimate`: The transformed estimate values
- `processinformation`: Information about the processes
- `estimationinformation`: Information about the estimation procedure

# Examples
```julia
# Take magnitude of a coherence estimate
mag_coh = abs(coherence_estimate)

# Extract real part of cross-spectrum
real_spec = real(cross_spectrum)

# Compute log-magnitude spectrum
log_spec = log(abs(spectrum))
```
"""
struct MarginallyTransformedEstimate{E, D, P, Q, N, S, A, T, IP, IE, F} <:
       AbstractEstimate{E, D, P, Q, N}
    argument::A
    estimate::T
    processinformation::IP
    estimationinformation::IE
    function MarginallyTransformedEstimate{E, S, N, F}(argument::A, estimate::T,
            processinfo::ProcessInformation{D}, estimationinfo::IE) where {
            E, S, N, F, A, T, D, IE}
        P, Q = checkinputs(argument, estimate, processinfo)
        IP = typeof(processinfo)
        return new{E, D, P, Q, N, S, A, T, IP, IE, F}(
            argument, estimate, processinfo, estimationinfo)
    end
end

function getestimatename(T::Type{<:MarginallyTransformedEstimate})::String
    return "$(gettransformname(T))($(getestimatename(getoriginaltype(T))))"
end
function getestimate(est::MarginallyTransformedEstimate)
    return est.estimate
end
function getargument(est::MarginallyTransformedEstimate)
    return est.argument
end

"""
    _construct_estimate_subset(::Type{<:MarginallyTransformedEstimate{...}}, trait, args...)

Internal constructor for creating estimate subsets with different traits.
"""
function _construct_estimate_subset(
        ::Type{<:MarginallyTransformedEstimate{E, D, P, Q, N, S, A, T, IP, IE, F}},
        trait::Type{<:EstimateTrait},
        args...
)::MarginallyTransformedEstimate where {E, D, P, Q, N, S, A, T, IP, IE, F}
    return MarginallyTransformedEstimate{trait, S, N, F}(args...)
end

# Type introspection utilities

"""
    getoriginaltype(::Type{<:MarginallyTransformedEstimate{E, D, P, Q, N, S}}) where {E, D, P, Q, N, S}

Extract the original estimate type before transformation.
"""
function getoriginaltype(::Type{<:MarginallyTransformedEstimate{
        E, D, P, Q, N, S}}) where {E, D, P, Q, N, S}
    return S
end

"""
    gettransformtype(est::MarginallyTransformedEstimate)

Get the type of transformation function applied to the estimate.
"""
gettransformtype(est::MarginallyTransformedEstimate) = gettransformtype(typeof(est))

function gettransformtype(::Type{<:MarginallyTransformedEstimate{
        E, D, P, Q, N, S, A, T, IP, IE, F}}) where {E, D, P, Q, N, S, A, T, IP, IE, F}
    return F
end

"""
    gettransformname(est::MarginallyTransformedEstimate)

Get a string representation of the transformation function name.
"""
gettransformname(est::MarginallyTransformedEstimate) = gettransformname(typeof(est))

function gettransformname(T::Type{<:MarginallyTransformedEstimate})::String
    return lstrip(string(nameof(gettransformtype(T))), '#')
end

# Core transformation functionality

"""
    apply_marginal_transform(transform::F, est::AbstractEstimate{E, D, P, Q, N}) where {F, E, D, P, Q, N}

Apply a transformation function element-wise to an estimate.

This is the core function for creating transformed estimates. It applies the given
transformation function to each element of the estimate while preserving all
metadata and argument structure.

# Arguments
- `transform::F`: A function to apply element-wise (e.g., `abs`, `real`, `log`)
- `est::AbstractEstimate`: The original estimate to transform

# Returns
A `MarginallyTransformedEstimate` wrapping the transformed values.

# Examples
```julia
# Apply absolute value transformation
abs_est = apply_marginal_transform(abs, complex_estimate)

# Apply logarithm with custom function
log_est = apply_marginal_transform(x -> log(x + 1e-10), spectrum)
```
"""
function apply_marginal_transform(
        transform::F, est::AbstractEstimate{E, D, P, Q, N}
)::MarginallyTransformedEstimate{E, D, P, Q, N} where {F, E, D, P, Q, N}
    argument = getargument(est)
    estimate = apply_marginal_transform(transform, getestimate(est))
    processinfo = getprocessinformation(est)
    estimationinfo = getestimationinformation(est)
    S = typeof(est)
    return MarginallyTransformedEstimate{E, S, N, F}(
        argument, estimate, processinfo, estimationinfo)
end

"""
    apply_marginal_transform(transform::F, x::AbstractArray{<:SMatrix}) where {F}

Apply transformation to arrays of static matrices.

For arrays containing static matrices (common in multi-process spectral estimates),
applies the transformation element-wise to each matrix element.
"""
function apply_marginal_transform(transform::F, x::AbstractArray{<:SMatrix}) where {F}
    return map(y -> transform.(y), x)
end

"""
    apply_marginal_transform(transform, x::AbstractArray{<:Number})

Apply transformation to arrays of numbers.

For arrays of scalar values, applies the transformation element-wise across the array.
"""
function apply_marginal_transform(transform, x::AbstractArray{<:Number})
    return transform.(x)
end

# Standard mathematical transformations
Base.real(x::AbstractEstimate) = apply_marginal_transform(real, x)
Base.imag(x::AbstractEstimate) = apply_marginal_transform(imag, x)
Base.conj(x::AbstractEstimate) = apply_marginal_transform(conj, x)
Base.abs(x::AbstractEstimate) = apply_marginal_transform(abs, x)
Base.abs2(x::AbstractEstimate) = apply_marginal_transform(abs2, x)
Base.angle(x::AbstractEstimate) = apply_marginal_transform(angle, x)
Base.log(x::AbstractEstimate) = apply_marginal_transform(log, x)
Base.exp(x::AbstractEstimate) = apply_marginal_transform(exp, x)
