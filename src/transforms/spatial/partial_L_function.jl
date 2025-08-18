struct PartialLFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    partial_L_function::T
    function PartialLFunction(radii::R, L::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, L)
        new{R,T,D,P}(radii, L)
    end
end

getargument(f::PartialLFunction) = f.radii
getestimate(f::PartialLFunction) = f.partial_L_function
getextrafields(::PartialLFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function partial_L_function(k::PartialKFunction{R,T,D,P}) where {R,T,D,P}
    @assert all(x -> all(y -> y ≥ 0, x), values(k.partial_K_function)) "Partial K function takes some negative values, so the L function cannot be used."
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    L = Dict(index => (val ./ V) .^ (1 / D) for (index, val) in k.partial_K_function)
    return PartialLFunction(k.radii, L, Val{D}())
end

function partial_L_function(k::PartialKFunction{R,T,2,P}) where {R,T,P}
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