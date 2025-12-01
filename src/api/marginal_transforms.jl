"""
    apply_marginal_transform(transform::F, est::AbstractEstimate{E, D, N}) where {F, E, D, N}

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
        transform::F, est::AbstractEstimate{E, D, N}; kwargs...) where {F, E, D, N}
    compute(MarginallyTransformedEstimate{E, D, typeof(est), F}, est; kwargs...)
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
