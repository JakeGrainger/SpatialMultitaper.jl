"""
	sdf2K(freq, spectra::AbstractArray{T,D}, r) where {D,T}

Computes the K function from the spectral density function.
"""
function sdf2K(freq, spectra::AbstractArray{T,D}, r, ::Nothing) where {D,T}
    prod(step, freq) * real(
        sum(
            f * sphere_weight(r, k, Val{D}()) for
            (f, k) in zip(spectra, Iterators.product(freq...))
        ),
    )
end

function sdf2K(freq, spectra::AbstractArray{T,D}, r, upscalefactor) where {T,D}
    itp = interpolate(spectra, BSpline(Linear()))
    upspec = itp(range.(1, size(spectra), upscalefactor .* (size(spectra) .- 1) .+ 1)...)
    upfreq = range.(first.(freq), last.(freq), size(upspec))
    sdf2K(upfreq, upspec, r, nothing)
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

"""
	partial_K(data, region; radii, nfreq, fmax, tapers, indices, upscalefactor)

Computes the partial K function from the `data` at radii `radii`.
Default is to compute this for all pairs of indices conditional on any index not included.
Alternatively, pass a vector of indices. If this is a vector of `Tuple{Int,Int}`, then this is computed partial on every other index.
If this is a `Tuple{Int,Int,AbstractVector{Int},AbstractVector{Int}}`, then this is computed partial on the specified indices, i.e.
The residual of `index[1]` partial `index[3]` with `index[2]` partial `index[4]`.
Setting `upscalefactor` to `nothing` does no upscaling. Otherwise, set to an `Int` or an `NTuple{D,Int}` to upscale the spectra for the computation of the partial K.
"""
function partial_K(
    data,
    region;
    radii,
    nfreq,
    fmax,
    tapers,
    indices = [(i, j) for i in eachindex(data), j in eachindex(data) if i <= j],
    upscalefactor = nothing,
)
    check_K_indices(indices, data)
    fhat = multitaper_estimate(data, region; tapers = tapers, nfreq = nfreq, fmax = fmax)
    zero_atom = atom_estimate.(data, Ref(region))
    return partial_K(fhat, zero_atom, radii, indices, upscalefactor)
end

function check_K_indices(indices, data)
    @assert indices isa
            AbstractVector{<:Tuple{Int,Int,AbstractVector{Int},AbstractVector{Int}}} ||
            indices isa AbstractVector{Tuple{Int,Int}}
    for index in indices
        for i in eachindex(index)
            @assert index[i] ⊆ eachindex(data) "problem with the indices, should all be subsets of $(eachindex(data)), but got $(index)"
        end
    end
end

function partial_K(
    fhat::SpectralEstimate,
    zero_atom,
    radii,
    indices::AbstractVector{Tuple{Int,Int}},
    upscalefactor = nothing,
)
    partial = partial_spectra(fhat)
    K = Dict(
        index => [
            sdf2K(
                partial.freq,
                partial[index[1], index[2]].partial_spectra .-
                (index[1] == index[2]) * zero_atom[index[1]],
                r,
                upscalefactor,
            ) for r in radii
        ] for index in indices
    )
    return (radii = radii, partial_K = K)
end

function partial_K(
    fhat::SpectralEstimate,
    zero_atom,
    radii,
    indices::AbstractVector{<:Tuple{Int,Int,AbstractVector{Int},AbstractVector{Int}}},
    upscalefactor = nothing,
)
    K = Dict(
        index =>
            [sdf2partialK(fhat, zero_atom, r, index, upscalefactor) for r in radii] for
        index in indices
    )
    return (radii = radii, partial_K = K)
end

function sdf2partialK(fhat, zero_atom, r, index, upscalefactor = nothing)
    pspec = partial_spectra(fhat, index[1], index[2], index[3], index[4])
    sdf2K(
        pspec.freq,
        pspec.partial_spectra .- (index[1] == index[2]) * zero_atom[index[1]],
        r,
        upscalefactor,
    )
end

##
atom_estimate(data::PointSet, region) = length(data) / unitless_measure(region)
atom_estimate(data::GeoTable, region) = atom_estimate(domain(data), values(data)[1], region)
atom_estimate(domain::CartesianGrid, rf, region) = 0.0
atom_estimate(domain::PointSet, marks, region) =
    error("Atom for the marked case is not yet implemented")
