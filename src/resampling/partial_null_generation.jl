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
    spec = multitaper_estimate(data, region; tapers = tapers, nfreq = nfreq, fmax = fmax)
    intensity = mean_estimate(data, region, mean_method)
    kernels = prediction_kernel(spec, radii = radii)
    kernel_integral =
        integrate_prediction_kernel.(Ref(radii), kernels.kernels, Val{embeddim(region)}())
    kernels_interp = ntuple(
        j -> ntuple(
            p -> linear_interpolation(
                radii,
                getindex.(kernels.kernels[j], p),
                extrapolation_bc = 0.0,
            ),
            Val{P - 1}(),
        ),
        Val{P}(),
    )
    ntuple(
        j -> create_single_intensity(
            j,
            intensity,
            kernel_integral[j],
            kernels_interp[j],
            data,
            grid,
        ),
        Val{P}(),
    )
end

function create_single_intensity(
    idx,
    intensities,
    kernel_integral,
    kernels_interp,
    data,
    grid,
)
    additional_processes = Not(idx)
    base_intensity =
        intensities[idx] -
        sum(kernel_integral .* collect(intensities)[additional_processes])
    sides = grid2side(grid)
    data_dep = collect(data)[additional_processes]
    intensity = [
        base_intensity + sum(
            sum(kernels_interp[j](norm(s .- unitless_coords(x))) for x in data_dep[j])
            for j in eachindex(data_dep)
        ) for s in Iterators.ProductIterator(sides)
    ]
    return georef((intensity = vec(intensity),), grid)
end

function integrate_prediction_kernel(radii, kernel, ::Val{D}) where {D}
    A = 2 * pi^(D / 2) / gamma(D / 2)
    return A * step(radii) * sum(k * r^(D - 1) for (k, r) in zip(kernel, radii))
end

function prediction_kernel(spec::SpectralEstimate; radii)
    kernel_ft = prediction_kernel_ft(spec)
    kernels = _ft2kernel.(Ref(kernel_ft.freq), kernel_ft.kernels, Ref(radii))
    return (radii = radii, kernels = kernels)
end

function _ft2kernel(freq::NTuple{D}, kernel_ft, radii) where {D}
    smooth_width = 5 # TODO: shouldn't be hardcoded
    return movingaverage(
        [
            prod(step, freq) * real(
                sum(
                    f * pcf_weight(radius, k, Val{D}()) for
                    (f, k) in zip(kernel_ft, Iterators.product(freq...))
                ),
            ) for radius in radii
        ],
        smooth_width,
    )
end

function movingaverage(x, width)
    y = similar(x)
    for i in eachindex(x)
        lo = max(1, i - div(width, 2))
        hi = min(length(x), i + div(width, 2))
        y[i] = zero(eltype(x))
        for idx = lo:hi
            y[i] += x[idx]
        end
        y[i] /= length(lo:hi)
    end
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