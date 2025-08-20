struct PartialCFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    partial_C_function::T
    function PartialCFunction(radii::R, C::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, C)
        new{R,T,D,P}(radii, C)
    end
end

getargument(f::PartialCFunction) = f.radii
getestimate(f::PartialCFunction) = f.partial_C_function
getextrafields(::PartialCFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function partial_C_function(
    f::PartialSpectra{D,F,P,N},
    zero_atom,
    radii,
) where {D,F,P,N}
    C = sdf2C(f, zero_atom, radii)
    return PartialCFunction(radii, C, Val{D}())
end

function partial_C_function(
    f::SpectralEstimate,
    zero_atom,
    radii,
    partial_type = UsualPartial()
)
    partial_C_function(partial_spectra(f, partial_type), zero_atom, radii)
end

function partial_C_function(
    data::Union{NTuple{P,Union{GeoTable,PointSet}},GeoTable,PointSet},
    region,
    radii;
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
    partial_type::PartialType = UsualPartial()
) where {P}
    fhat = multitaper_estimate(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        mean_method = mean_method,
    )
    zero_atom = atom_estimate(data, region)
    return partial_C_function(fhat, zero_atom, radii, partial_type)
end