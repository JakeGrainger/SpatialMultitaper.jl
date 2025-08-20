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
    indices = default_indices(f),
) where {D,F,P,N}
    check_spectral_indices(indices, f)
    C = sdf2C(f, zero_atom, radii, indices)
    return PartialCFunction(radii, C, Val{D}())
end

function partial_C_function(
    f::SpectralEstimate,
    zero_atom,
    radii,
    indices = default_indices(f),
)
    check_spectral_indices(indices, f)
    partial_C_function(partial_spectra(f), zero_atom, radii, indices)
end

function partial_C_function(
    f::SpectralEstimate{D,F,P,N},
    zero_atom,
    radii,
    indices::AbstractVector{<:Tuple{Int,Int,AbstractVector{Int},AbstractVector{Int}}},
) where {D,F,P,N}
    check_spectral_indices(indices, f)
    C = Dict(
        index => sdf2C(
            partial_spectra(f, index[1], index[2], index[3], index[4]),
            (index[1] == index[2]) * zero_atom[index[1]],
            radii,
        ) for index in indices
    )
    return PartialCFunction(radii, C, Val{D}())
end

function partial_C_function(
    data::Union{NTuple{P,Union{GeoTable,PointSet}},GeoTable,PointSet},
    region,
    radii,
    indices = default_indices(data);
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
) where {P}
    check_spectral_indices(indices, data)
    fhat = multitaper_estimate(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        mean_method = mean_method,
    )
    zero_atom = atom_estimate(data, region)
    return partial_C_function(fhat, zero_atom, radii, indices)
end