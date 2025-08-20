struct LFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    L_function::T
    function LFunction(radii::R, L::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, L)
        new{R,T,D,P}(radii, L)
    end
end

getargument(f::LFunction) = f.radii
getestimate(f::LFunction) = f.L_function
getextrafields(::LFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function K2L(k_function, ::Val{D}) where {D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    return [sign.(k) .* (abs.(k) ./ V) .^ (1 / D) for k in k_function] # broadcast because may contain SVectors
end

function K2L(k_function, ::Val{2})
    return [sign.(k) .* sqrt.(abs.(k) ./ pi) for k in k_function] # broadcast because may contain SVectors
end

function L_function(k::KFunction{R,T,D,P}) where {R,T,D,P}
    return LFunction(k.radii, K2L(k.K_function, Val{D}()), Val{D}())
end

function L_function(
    data,
    region,
    radii;
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    k = K_function(
        data,
        region,
        radii;
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
    )
    return L_function(k)
end