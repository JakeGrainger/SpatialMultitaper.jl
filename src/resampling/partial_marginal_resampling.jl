struct PartialMarginalResampler{T, G}
    Λ::T
    region::G
end

Base.getindex(m::PartialMarginalResampler, i) = PartialMarginalResampler(m.Λ[i], m.region)

function PartialMarginalResampler(
        data::NTuple{P, PointSet},
        region;
        tapers,
        nfreq,
        fmax
) where {P}
    Λ = create_intensities(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax
    )
    return PartialMarginalResampler(Λ, region)
end

## sampling
function Base.rand(
        rng::AbstractRNG, m::PartialMarginalResampler{NTuple{P, T}}) where {P, T}
    ntuple(p -> Base.rand(rng, m[p]), Val{P}())
end

function Base.rand(rng::AbstractRNG, m::PartialMarginalResampler{<:GeoTable})
    Λ, region = m.Λ, m.region
    @assert Λ isa GeoTable "you may need to index your PartialMarginalResampler to get a single intensity, currently it is likely a collection of $(length(Λ)) intensities"
    λ = reshape(first(values(Λ)), size(domain(Λ)))

    λ₀ = maximum(l for l in λ if !isnan(l))
    grid_min = unitless_minimum(domain(Λ))
    grid_spacing = unitless_spacing(domain(Λ))
    proposal = rand(rng, PoissonProcess(λ₀), region)
    thinned = eltype(proposal)[]
    for p in proposal
        pcoords = SpatialMultitaper.unitless_coords(p)
        if rand(rng) ≤ λ[intensity_index(pcoords, grid_min, grid_spacing)] / λ₀
            push!(thinned, p)
        end
    end
    mask(Meshes.PointSet(thinned), region)
end

function intensity_index(pcoords, grid_min, grid_spacing)
    CartesianIndex(
        @. ceil(Int, (pcoords - grid_min) / grid_spacing)
    )
end

## construction
function create_intensities(
        data::NTuple{P, PointSet},
        region;
        tapers,
        nfreq,
        fmax,
        mean_method::MeanEstimationMethod = DefaultMean()
) where {P}
    spec = multitaper_estimate(data, region; tapers = tapers, nfreq = nfreq, fmax = fmax)
    intensity = mean_estimate(data, region, mean_method)
    data_ft_full = fft_only.(data, Ref(region), nfreq = nfreq, fmax = fmax)
    data_ft = zeros(SVector{P, eltype(first(data_ft_full))}, size(first(data_ft_full)))
    for i in eachindex(first(data_ft_full))
        data_ft[i] = SVector{P, eltype(first(data_ft_full))}(getindex.(
            data_ft_full, Ref(i)))
    end

    kernel_ft = prediction_kernel_ft(spec)
    ntuple(
        j -> create_single_intensity(
            j,
            intensity,
            data_ft,
            kernel_ft
        ),
        Val{P}()
    )
end

function create_single_intensity(
        idx, intensities::SVector{P, T}, data_ft, kernel_ft) where {P, T}
    freq = kernel_ft.freq
    ker = kernel_ft.kernels[idx]
    λz = SVector{P - 1, T}((intensities[1:(idx - 1)]..., intensities[(idx + 1):P]...))

    base = intensities[idx] - (ker[findfirst.(Ref(iszero), freq)...] * λz)[1]

    idx_other = static_not(Val{P}(), idx)
    adjustment_ft = getindex.(ker .* getindex.(data_ft, Ref(idx_other)), 1)
    adjustment = bfft(adjustment_ft)
    adjustment .*= prod(step, freq)

    intensity = abs.(base .+ (adjustment)) # TODO: abs is probably not ideal
    grid = CartesianGrid(
        length.(freq), ntuple(zero, length(freq)), 1 ./
                                                   (length.(freq) .* step.(freq)))
    georef((intensity = vec(intensity),), grid)
end

"""
    fft_only(points::PointSet, region; nfreq, fmax)

Does a non-uniform FFT of points, but moves them so the origin of the bounding box of the region is at zero.
Makes no adjustments for tapering etc.
"""
function fft_only(points::PointSet, region; nfreq, fmax)
    # bbox = boundingbox(region)
    # translation = Translate(.-Meshes.to(minimum(bbox)).coords)
    # t_points = translation(points) # translate to origin
    # t_region = translation(region)
    # shifts will cancel out due to the other transforms also having the same shift
    out = nufft_anydomain(
        region, nfreq, fmax, points, ones(ComplexF64, length(points)), -1, 1e-14) # TODO: make work for grid data
    return reshape(out, size(out)[1:(end - 1)]) # last dimension is singletons
end

"""
    prediction_kernel_ft(spec)

Computes the Fourier transform of the prediction kernel given estimates of the spectral density function.
"""
function prediction_kernel_ft(spec::SpectralEstimate{
        F, N, I, T, D, P, P}) where {F, N, I, T, D, P}
    kernels = ntuple(
        idx -> apply_transform(single_prediction_kernel_ft, spec.power,
            idx, getestimationinformation(spec).ntapers),
        Val{P}()
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
        S_mat::SMatrix{P, P, T, L},
        idx::Int,
        ::Nothing
) where {P, T, L}
    @boundscheck checkbounds(1:P, idx)
    other_idx = StaticArrays.sacollect(
        SVector{P - 1, Int}, ApplyArray(vcat, 1:(idx - 1), (idx + 1):P))
    S_XZ = S_mat[SVector{1, Int}(idx), other_idx]
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
