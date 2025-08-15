struct LFunction{R,T,D}
    radii::R
    L_function::T
    LFunction(radii::R, L::T, ::Val{D}) where {R,T,D} = new{R,T,D}(radii, L)
end

function L_function(k::KFunction{R,T,D}) where {R,T,D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    L = Dict(index => (val ./ V) .^ (1 / D) for (index, val) in k.K_function)
    return LFunction(k.radii, L, Val{D}())
end

function L_function(k::KFunction{R,T,2}) where {R,T}
    L = Dict(index => sqrt.(val ./ pi) for (index, val) in k.K_function)
    return LFunction(k.radii, L, Val{2}())
end

function L_function(
    data,
    region,
    radii,
    indices = default_indices(data);
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    k = K_function(
        data,
        region,
        radii,
        indices;
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
    )
    return L_function(k)
end