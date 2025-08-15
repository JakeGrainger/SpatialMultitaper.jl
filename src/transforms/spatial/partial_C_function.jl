struct PartialCFunction{R,T,D}
    radii::R
    C::T
    PartialCFunction(radii::R, C::T, ::Val{D}) where {R,T,D} = new{R,T,D}(radii, C)
end

function partial_C_function(
    f::PartialSpectra{D,F,P,N},
    zero_atom,
    radii,
    indices::AbstractVector{Tuple{Int,Int}} = default_indices(f),
) where {D,F,P,N}
    check_spectral_indices(indices, f)
    C = Dict(
        index => [
            sdf2C(f[index...], zero_atom[index[1]] * (index[1] == index[2]), r) for
            r in radii
        ] for index in indices
    )
    return PartialCFunction(radii, C, Val{D}())
end

function partial_C_function(
    f::SpectralEstimate,
    zero_atom,
    radii,
    indices::AbstractVector{Tuple{Int,Int}} = default_indices(f),
)
    check_spectral_indices(indices, f)
    partial_C_function(partial_spectra(f), zero_atom, radii; indices = indices)
end

function partial_C_function(
    f::SpectralEstimate{D,F,P,N},
    zero_atom,
    radii,
    indices::AbstractVector{<:Tuple{Int,Int,AbstractVector{Int},AbstractVector{Int}}},
) where {D,F,P,N}
    check_spectral_indices(indices, f)
    C = Dict(
        index => [
            sdf2K(
                partial_spectra(f, index[1], index[2], index[3], index[4]),
                (index[1] == index[2]) * zero_atom[index[1]],
                r,
            ) for r in radii
        ] for index in indices
    )
    return PartialCFunction(radii, C, Val{D}())
end

function partial_C_function(
    data,
    region,
    radii,
    indices::AbstractVector{Tuple{Int,Int}} = default_indices(data);
    nfreq,
    fmax,
    tapers,
)
    check_spectral_indices(indices, data)
    fhat = multitaper_estimate(data, region; tapers = tapers, nfreq = nfreq, fmax = fmax)
    zero_atom = atom_estimate.(data, Ref(region))
    return partial_C_function(fhat, zero_atom, radii; indices = indices)
end
