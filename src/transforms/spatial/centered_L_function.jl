struct CenteredLFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    centered_L_function::T
    function CenteredLFunction(radii::R, L::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, L)
        new{R,T,D,P}(radii, L)
    end
end

getargument(f::CenteredLFunction) = f.radii
getestimate(f::CenteredLFunction) = f.centered_L_function
getextrafields(::CenteredLFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function L2centeredL(radii, k)
    k .- radii
end

function centered_L_function(l::LFunction{R,T,D,1}) where {R,T,D}
    return CenteredLFunction(l.radii, L2centeredL(l.radii, l.L_function), Val{D}())
end

function centered_L_function(k::LFunction{R,T,D,P}) where {R,T,D,P}
    L = Dict(index => L2centeredL(k.radii, val) for (index, val) in k.L_function)
    return CenteredLFunction(k.radii, L, Val{D}())
end

centered_L_function(k::KFunction) = centered_L_function(L_function(k))

function centered_L_function(
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
    return centered_L_function(k)
end