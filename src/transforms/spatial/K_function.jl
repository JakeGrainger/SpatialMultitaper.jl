struct KFunction{R,T,D}
    radii::R
    K::T
    KFunction(radii::R, K::T, ::Val{D}) where {R,T,D} = new{R,T,D}(radii, K)
end

function K_function(c::CFunction{R,T,D}, λ) where {R,T,D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    K = Dict(
        index => val ./ λ[index[1]] * λ[index[2]] .+ c.radii .^ d .* V for
        (index, val) in c.C
    )
    return KFunction(c.radii, K, Val{D}())
end

function K_function(
    data,
    region,
    indices = default_indices(data);
    radii,
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    c = C_function(
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
    return K_function(c, λ)
end