struct PartialKFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    partial_K_function::T
    function PartialKFunction(radii::R, K::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, K)
        new{R,T,D,P}(radii, K)
    end
end

getargument(f::PartialKFunction) = f.radii
getestimate(f::PartialKFunction) = f.partial_K_function
getextrafields(::PartialKFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function partial_K_function(c::PartialCFunction{R,T,D,1}, λ) where {R,T,D}
    invλ = 1 / λ[1]
    return PartialKFunction(
        c.radii,
        C2K(c.radii, c.partial_C_function, invλ, invλ, Val{D}()),
        Val{D}(),
    )
end

function partial_K_function(
    c::PartialCFunction{R,T,D,P},
    λ::NTuple{P,<:Number},
) where {R,T<:AbstractArray,D,P}
    invλ = diagm(1 ./ SVector(λ...))
    return PartialKFunction(
        c.radii,
        C2K(c.radii, c.partial_C_function, invλ, invλ, Val{D}()),
        Val{D}(),
    )
end

function partial_K_function(c::PartialCFunction{R,T,D,P}, λ) where {R,T<:Dict,D,P}
    K = Dict(
        index => C2K(c.radii, val, 1 / λ[index[1]], 1 / λ[index[2]], Val{D}()) for
        (index, val) in c.partial_C_function
    )
    return PartialKFunction(c.radii, K, Val{D}())
end

"""
	partial_K_function(data, region, radii, indices; nfreq, fmax, tapers, mean_method)

Computes the partial K function from the `data` at radii `radii`.
Default is to compute this for all pairs of indices conditional on any index not included.
Alternatively, pass a vector of indices. If this is a vector of `Tuple{Int,Int}`, then this is computed partial on every other index.
If this is a `Tuple{Int,Int,AbstractVector{Int},AbstractVector{Int}}`, then this is computed partial on the specified indices, i.e.
The residual of `index[1]` partial `index[3]` with `index[2]` partial `index[4]`.
"""
function partial_K_function(
    data,
    region,
    radii,
    indices = default_indices(data);
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    c = partial_C_function(
        data,
        region,
        radii,
        indices;
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
    )
    λ = mean_estimate(data, region, mean_method)
    return partial_K_function(c, λ)
end