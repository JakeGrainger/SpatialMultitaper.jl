struct MarginalResampler{T,G}
    Λ::T
    region::G
end

Base.getindex(m::MarginalResampler, i) = MarginalResampler(m.Λ[i], m.region)

function MarginalResampler(
    data::NTuple{P,PointSet},
    region;
    tapers,
    nfreq,
    fmax,
    radii,
    grid,
) where {P}
    Λ = create_intensities(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        radii = radii,
        grid = grid,
    )
    return MarginalResampler(Λ, region)
end

function Base.rand(m::MarginalResampler)
    Λ, region = m.Λ, m.region
    @assert Λ isa GeoTable "you may need to index your MarginalResampler to get a single intensity, currently it is likely a collection of $(length(Λ)) intensities"
    λ = reshape(first(values(Λ)), size(domain(Λ)))
    pts = grid2side(domain(Λ))

    λ₀ = maximum(l for l in λ if !isnan(l))
    grid_min = first.(pts)
    proposal = rand(PoissonProcess(λ₀), region)
    thinned = eltype(proposal)[]
    for p ∈ proposal
        pcoords = SpatialMultitaper.unitless_coords(p)
        intensity_index = CartesianIndex(
            @. min(floor(Int, (pcoords - grid_min) / step(pts)) + 1, length(pts))
        )
        if rand() ≤ λ[intensity_index] / λ₀
            push!(thinned, p)
        end
    end
    mask(Meshes.PointSet(thinned), region)
end

"""
    create_intensities(data::NTuple{P,PointSet}, region; tapers, nfreq, fmax, radii, grid) where {P}

Constructs the internal intensities for each of the processes partial on all the others.

``\\lambda_{X\\cdot Z}(u) = \\lambda_X - \\sum_{p} \\lambda_{Z_p} \\int \\psi_j(x)d x + \\sum_{p} \\sum_{x \\in Z_p} \\psi_j(u - x)``

where ``\\tilde\\psi_j(k) = [f_{XZ}(k)f_{ZZ}(k)^{-1}]_j`` is the Fourier transform of the jth prediction kernel.
"""
function create_intensities(
    data::NTuple{P,PointSet},
    region;
    tapers,
    nfreq,
    fmax,
    radii,
    grid,
    mean_method::MeanEstimationMethod = DefaultMean(),
) where {P}
    mask.(data, Ref(region))
    data, dim = check_spatial_data(data)
    mean_method = check_mean_method(mean_method, data)
    J_n = tapered_dft(data, tapers, nfreq, fmax, region, mean_method)
    freq = make_freq(nfreq, fmax, dim)
    power = dft2spectralmatrix(J_n)
    spec = SpectralEstimate(freq, power, length(tapers))
    intensities = mean_estimate.(data, Ref(region), mean_method)
    kernels_ft = prediction_kernel_ft(spec)
    return ntuple(
        idx -> create_single_intensity(idx, intensities, kernels_ft, J_n),
        Val{P}(),
    )
end

function create_single_intensity(idx, intensities, kernels_ft, J_n::NTuple{P}) where {P}
    freq = kernels_ft.freq
    λ = intensities[idx]
    other_idx =
        StaticArrays.sacollect(SVector{P - 1,Int}, ApplyArray(vcat, 1:idx-1, idx+1:P))
    freq_version = [
        sum(ψ[j] * J_n[other_idx[j]][i] for j = 1:(P-1)) for
        (i, ψ) in enumerate(kernels_ft.kernels[idx])
    ]
    intensity_partial =
        λ .+
        (length(freq_version) * prod(step.(freq))) .*
        fftshift(ifft(ifftshift(freq_version)))
    grid_sides = fftshift.(fftfreq.(length.(freq), length.(freq) ./ step.(freq)))
    return georef((intensity = vec(intensity_partial),), side2grid(grid_sides))
end

# commented out is for anisotropic version
# function prediction_kernel_ft2space(freq, power::AbstractArray{D,T}) where {D,T<:Number}
#     @assert length(freq) == ndims(power) - 1
#     (length(power) * prod(step.(freq))) .*
#     fftshift(ifft(ifftshift(power, 2:ndims(power)), 2:ndims(power)), 2:ndims(power))
# end

# function prediction_kernel_ft2space(freq, power::AbstractArray{S,SMatrix})
#     prediction_kernel_ft2space(freq, svectors2array(power))
# end

# function svectors2array(x::Array{D,SVector{P,T}}) where {D,P,T}
#     reinterpret(reshape, T, x)
# end


"""
    prediction_kernel_ft(spec)

Computes the Fourier transform of the prediction kernel given estimates of the spectral density function.
"""
function prediction_kernel_ft(spec::SpectralEstimate{D,F,P,N,T}) where {D,F,P,N,T}
    kernels = ntuple(
        idx ->
            apply_transform(single_prediction_kernel_ft, spec.power, idx, spec.ntapers),
        Val{P}(),
    )
    return (freq = spec.freq, kernels = kernels)
end

function single_prediction_kernel_ft(S_mat::AbstractMatrix, idx::Int, ntapers::Int)
    scaling = ntapers / (ntapers - size(S_mat, 1) + 1)
    return single_prediction_kernel_ft(S_mat, idx, nothing) .* scaling
end

"""
    single_prediction_kernel_ft(S_mat::SMatrix{P,P,T,L}, idx::Int, ::Nothing)

Computes the Fourier transform of the prediction kernel for a single variable given an estimate of the spectral density function.
"""
function single_prediction_kernel_ft(
    S_mat::SMatrix{P,P,T,L},
    idx::Int,
    ::Nothing,
) where {P,T,L}
    @boundscheck checkbounds(1:P, idx)
    other_idx =
        StaticArrays.sacollect(SVector{P - 1,Int}, ApplyArray(vcat, 1:idx-1, idx+1:P))
    S_XZ = S_mat[SVector{1,Int}(idx), other_idx]
    S_ZZ = S_mat[other_idx, other_idx]
    return S_XZ / S_ZZ
end

"""
    single_prediction_kernel_ft(S_mat::AbstractMatrix, idx::Int, ::Nothing)

Computes the Fourier transform of the prediction kernel for a single variable given an estimate of the spectral density function.
"""
function single_prediction_kernel_ft(S_mat::AbstractMatrix, idx::Int, ::Nothing)
    @boundscheck checkbounds(1:size(S_mat, 1), idx)
    S_XZ = S_mat[idx, Not(idx)]
    S_ZZ = S_mat[Not(idx), Not(idx)]
    return S_XZ / S_ZZ
end