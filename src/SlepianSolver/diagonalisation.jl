"""
	make_concentration_operator(R::AbstractArray,K::AbstractArray,real_tapers=Val{true}())

Constructs the concentration operator. Provide `real_tapers=Val{true}()` for real valued tapers.
Choose `real_tapers=Val{false}()` for complex valued tapers.
"""
function make_concentration_operator(
    R::AbstractArray,
    K::AbstractArray,
    real_tapers = Val{true}(),
)
    size(R) == size(K) || error("Regions should be the same size in space and freq.")
    P = SupportProjection(R)
    L = SupportProjection(K)
    Q = FourierTransform(size(inputspace(P)))
    if real_tapers isa Val{true}
        V = Vec(size(R), Val{Float64}())
        Re = RealPart{Float64}()
        operator = V * P * Re * inv(Q) * L * Q * P * inv(V)
    else
        V = Vec(size(R), Val{ComplexF64}())
        operator = V * P * inv(Q) * L * Q * P * inv(V)
    end
    return operator, inv(V)
end

"""
	compute_eigenfunctions(R::AbstractArray,K::AbstractArray,ntapers=1;imag_tol=1e-5,real_tapers=true, tol=0.0)

Computes the eigenfunctions of the concentration operator. Provide `real_tapers=true` for real valued tapers.
"""
function compute_eigenfunctions(
    R::AbstractArray{Bool,D},
    K::AbstractArray{Bool,D},
    ntapers = 1;
    imag_tol = 1e-5,
    real_tapers = true,
    tol = 0.0,
) where {D}
    _compute_eigenfunctions(R, K, ntapers, imag_tol, Val{real_tapers}(), tol)
end
function _compute_eigenfunctions(
    R::AbstractArray{Bool,D},
    K::AbstractArray{Bool,D},
    ntapers,
    imag_tol,
    real_tapers,
    tol,
) where {D}
    @assert size(K) == size(R)
    conc, reshp = make_concentration_operator(R, K, real_tapers)
    λ, vh = eigs(conc, nev = ntapers, which = :LR, tol = tol)
    maximum(imag.(vh)) < imag_tol || @warn "Imaginary part was larger than expected."
    if real_tapers isa Val{true}
        h = [reshp * real.(vh[:, i]) for i in axes(vh, 2)]
    else
        h = [reshp * (vh[:, i]) for i in axes(vh, 2)]
    end
    return h, λ
end
