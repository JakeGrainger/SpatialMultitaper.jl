"""
    sdf2K(freq, spectra::AbstractArray{T,D}, r) where {D,T}

Computes the K function from the spectral density function.
"""
function sdf2K(freq, spectra::AbstractArray{T,D}, r) where {D,T}
    prod(step, freq) * real(
        sum(
            f * sphere_weight(r, k, Val{D}()) for
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

function partial_K(data, R, tapers; region, nfreq, fmax, indices = [[i, j] for i in eachindex(data), j in eachindex(data) if i <= j])
    @assert all(length(index)==2 && 1 ≤ index[1] ≤ length(data) && 1 ≤ index[2] ≤ length(data) for index in indices)
    fhat = multitaper_estimate(data, tapers; region = region, nfreq = nfreq, fmax = fmax)
    zero_atom = atom_estimate.(data, Ref(region))
    return partial_K(fhat, zero_atom, R, indices)
end

function partial_K(fhat::SpectralEstimate, zero_atom, R, indices)
    partial = partial_spectra(fhat)
    return R, Dict(index => [
        sdf2K(partial.freq, partial.partial_spectra[index[1], index[2], :, :] .- (index[1]==index[2]) * zero_atom[index[1]], r) for r in R
    ] for index in indices)
end

##
atom_estimate(data::PointSet, region) = length(data)/unitless_measure(region)
atom_estimate(data::GeoTable, region) = atom_estimate(domain(data), values(data)[1], region)
atom_estimate(domain::CartesianGrid, rf, region) = 0.0
atom_estimate(domain::PointSet, marks, region) =  error("Atom for the marked case is not yet implemented")
