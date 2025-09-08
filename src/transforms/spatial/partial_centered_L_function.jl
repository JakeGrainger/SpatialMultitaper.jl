struct PartialCenteredLFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    partial_centered_L_function::T
    function PartialCenteredLFunction(radii::R, L::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, L)
        new{R,T,D,P}(radii, L)
    end
end

getargument(f::PartialCenteredLFunction) = f.radii
getestimate(f::PartialCenteredLFunction) = f.partial_centered_L_function
getextrafields(::PartialCenteredLFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function partial_centered_L_function(l::PartialLFunction{R,T,D,P}) where {R,T,D,P}
    return PartialCenteredLFunction(
        l.radii,
        L2centeredL(l.radii, l.partial_L_function),
        Val{D}(),
    )
end

partial_centered_L_function(k::PartialKFunction) =
    partial_centered_L_function(partial_L_function(k))

function partial_centered_L_function(
    data,
    region;
    radii,
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
    partial_type::PartialType = UsualPartial(),
    freq_radii = default_rotational_radii(nfreq, fmax),
    rotational_method = default_rotational_kernel(freq_radii),
)
    k = partial_K_function(
        data,
        region;
        radii = radii,
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
        partial_type = partial_type,
        freq_radii = freq_radii,
        rotational_method = rotational_method,
    )
    return partial_centered_L_function(k)
end