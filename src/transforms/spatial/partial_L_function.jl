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
    return PartialLFunction(k.radii, K2L(k.partial_K_function, Val{D}()), Val{D}())
end

function partial_L_function(
    data,
    region,
    radii;
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
    partial_type::PartialType = UsualPartial()
)
    k = partial_K_function(
        data,
        region,
        radii;
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
        partial_type = partial_type
    )
    return partial_L_function(k)
end