struct PartialKFunction{R,T,D}
    radii::R
    K::T
    PartialKFunction(radii::R, K::T, ::Val{D}) where {R,T,D} = new{R,T,D}(radii, K)
end

function partial_K_function(c::PartialCFunction{R,T,D}, λ) where {R,T,D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    K = Dict(
        index => val ./ λ[index[1]] * λ[index[2]] .+ c.radii .^ d .* V for
        (index, val) in c.C
    )
    return PartialKFunction(c.radii, K, Val{D}())
end

"""
	partial_K(data, region, radii, indices; nfreq, fmax, tapers, mean_method)

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
        indices;
        radii = radii,
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
    )
    λ = mean_estimate(data, region, mean_method)
    return partial_K_function(c, λ)
end