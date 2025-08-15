struct CFunction{R,T,D}
    radii::R
    C::T
    CFunction(radii::R, C::T, ::Val{D}) where {R,T,D} = new{R,T,D}(radii, C)
end

function C_function(
    f::SpectralEstimate{D,F,P,N},
    zero_atom,
    radii,
    indices = default_indices(f),
) where {D,F,P,N}
    check_spectral_indices(indices, f)
    C = Dict(
        index => [
            sdf2C(f[index...], zero_atom[index[1]] * (index[1] == index[2]), r) for
            r in radii
        ] for index in indices
    )
    return CFunction(radii, C, Val{D}())
end

function C_function(
    data,
    region,
    radii,
    indices = default_indices(data);
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    check_spectral_indices(indices, data)
    fhat = multitaper_estimate(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        mean_method = mean_method,
    )
    zero_atom = atom_estimate.(data, Ref(region))
    return C_function(fhat, zero_atom, radii; indices = indices)
end

"""
    sdf2C(f, zero_atom, radius)

Takes some form of spectra and returns the C function for the `radius`.
"""
function sdf2C(
    f::Union{SpectralEstimate{D,F,P,N},PartialSpectra{D,F,P,N}},
    zero_atom,
    radius::Number,
) where {D,F,P,N}
    freq = getfreq(f)
    spectra = getestimate(f) .- zero_atom
    prod(step, freq) * real(
        sum(
            f * sphere_weight(radius, k, Val{D}()) for
            (f, k) in zip(spectra, Iterators.product(freq...))
        ),
    )
end

function sphere_weight(r, u, ::Val{D}) where {D}
    x = norm(u)
    if D === 1
        return 2r * sinc(2r * x)
    elseif D === 2
        return x < 1e-10 ? pi * r^2 : (r / x) * besselj1(2π * r * x)
    else
        return (r / x)^(D / 2) * besselj(D / 2, 2π * r * x)
    end
end