using LinearAlgebra
using BandedMatrices
using BSplineKit

import BSplineKit: Natural

"""
    spar_to_lambda(order::BSplineOrder, xs::AbstractVector, ys::AbstractVector, spar::Real;
                   weights::Union{Nothing, AbstractVector} = nothing)

Convert the smoothing parameter `spar` to the smoothing parameter `λ` for a cubic B-spline

See R's `smooth.spline` documentation for details on the conversion.
```
"""
function spar_to_lambda(
        order::BSplineOrder{4}, xs::AbstractVector, ys::AbstractVector, spar::Real;
        weights::Union{Nothing, AbstractVector} = nothing)
    ## start code from BPlineKit (should be pull request in that package)
    eachindex(xs) == eachindex(ys) ||
        throw(DimensionMismatch("x and y vectors must have the same length"))
    N = length(xs)
    cs = similar(ys)

    T = eltype(xs)

    # Create natural cubic B-spline basis with knots = input points
    B = BSplineBasis(order, copy(xs))
    R = RecombinedBSplineBasis(B, Natural())

    # Compute collocation matrices for derivatives 0 and 2
    A = BandedMatrix{T}(undef, (N, N), (1, 1))  # 3 bands are enough
    D = similar(A)
    collocation_matrix!(A, R, xs, Derivative(0))
    collocation_matrix!(D, R, xs, Derivative(2))

    # Matrix containing grid steps (banded + symmetric, so we don't compute the lower part)
    Δ_upper = BandedMatrix{T}(undef, (N, N), (0, 1))
    fill!(Δ_upper, 0)
    @inbounds for i in axes(Δ_upper, 1)[2:(end - 1)]
        # Δ_upper[i, i - 1] = xs[i] - xs[i - 1]  # this one is obtained by symmetry
        Δ_upper[i, i] = 2 * (xs[i + 1] - xs[i - 1])
        Δ_upper[i, i + 1] = xs[i + 1] - xs[i]
    end
    Δ = Hermitian(Δ_upper)  # symmetric matrix with 3 bands

    # The integral of the squared second derivative is (H * cs) ⋅ (D * cs) / 6.
    H = Δ * D

    # Directly construct LHS matrix
    # M = Hermitian((H'D + D'H) * (λ / 6) + 2 * A' * W * A)  # usually positive definite

    # Construct LHS matrix trying to reduce computations
    B = H' * D  # this matrix has 5 bands
    ## finish code from BPlineKit

    # Compute tr(Σ) = tr(H'D + D'H) / 6
    # Since the integral of squared second derivative is (H * cs) ⋅ (D * cs) / 6
    if weights === nothing
        tr_XWX = tr(A' * A)
    else
        tr_XWX = tr(A' * W * A)
    end

    penalty_mat = (B .+ B') ./ 6
    tr_Sigma = tr(penalty_mat)

    # Compute ratio
    ratio = tr_XWX / tr_Sigma

    # Convert spar to lambda
    return ratio * 256.0^(3.0 * spar - 1)
end
