struct PartialLFunction{R,T,D}
    radii::R
    partial_L_function::T
    PartialLFunction(radii::R, L::T, ::Val{D}) where {R,T,D} = new{R,T,D}(radii, L)
end

function partial_L_function(k::PartialKFunction{R,T,D}) where {R,T,D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    L = Dict(index => (val ./ V) .^ (1 / D) for (index, val) in k.partial_K_function)
    return PartialLFunction(k.radii, L, Val{D}())
end

function partial_L_function(k::PartialKFunction{R,T,2}) where {R,T}
    L = Dict(index => sqrt.(val ./ pi) for (index, val) in k.partial_K_function)
    return PartialLFunction(k.radii, L, Val{2}())
end

function partial_L_function(
    data,
    region,
    radii,
    indices = default_indices(data);
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    k = partial_K_function(
        data,
        region,
        radii,
        indices;
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
    )
    return partial_L_function(k)
end