struct PartialCenteredLFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    partial_partial_centered_L_function::T
    function PartialCenteredLFunction(radii::R, L::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, L)
        new{R,T,D,P}(radii, L)
    end
end

getargument(f::PartialCenteredLFunction) = f.radii
getestimate(f::PartialCenteredLFunction) = f.partial_centered_L_function
getextrafields(::PartialCenteredLFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function partial_centered_L_function(l::PartialLFunction{R,T,D,P}) where {R,T,D,P}
    return PartialCenteredLFunction(l.radii, L2centeredL(l.radii, l.L_function), Val{D}())
end

function partial_centered_L_function(k::PartialLFunction{R,T,D,P}) where {R,T<:Dict,D,P}
    L = Dict(index => L2centeredL(k.radii, val) for (index, val) in k.L_function)
    return PartialCenteredLFunction(k.radii, L, Val{D}())
end

partial_centered_L_function(k::PartialKFunction) =
    partial_centered_L_function(L_function(k))

function partial_centered_L_function(
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
    return partial_centered_L_function(k)
end