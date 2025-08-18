struct KFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    K_function::T
    function KFunction(radii::R, K::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, K)
        new{R,T,D,P}(radii, K)
    end
end

getargument(f::KFunction) = f.radii
getestimate(f::KFunction) = f.K_function
getextrafields(::KFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function C2K(radii, c, λ1, λ2, ::Val{D}) where {D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    return c ./ (λ1 * λ2) .+ radii .^ D .* V
end


function K_function(c::CFunction{R,T,D,1}, λ) where {R,T,D}
    return KFunction(c.radii, C2K(c.radii, c.C_function, λ[1], λ[1], Val{D}()), Val{D}())
end

function K_function(c::CFunction{R,T,D,P}, λ) where {R,T,D,P}
    K = Dict(
        index => C2K(c.radii, val, λ[index[1]], λ[index[2]], Val{D}()) for
        (index, val) in c.C_function
    )
    return KFunction(c.radii, K, Val{D}())
end

function K_function(
    data,
    region,
    radii,
    indices = default_indices(data);
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    c = C_function(
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
    return K_function(c, λ)
end