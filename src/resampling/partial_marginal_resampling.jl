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
        fmax,
        radii,
        grid,
        smooth_width = nothing
) where {P}
    Λ = create_intensities(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        radii = radii,
        grid = grid,
        smooth_width = smooth_width
    )
    return PartialMarginalResampler(Λ, region)
end

## sampling
function Base.rand(rng::AbstractRNG, m::PartialMarginalResampler)
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
        idx, intensities::NTuple{P, T}, data_ft, kernel_ft) where {P, T}
    freq = kernel_ft.freq
    ker = kernel_ft.kernels[idx]
    λz = SVector{P - 1, T}((intensities[1:(idx - 1)]..., intensities[(idx + 1):P]...))

    base = intensities[idx] - ker[findfirst.(Ref(iszero), freq)] * λz

    idx_other = static_not(Val{P}(), idx)
    adjustment_ft = kernel_ft .* getindex.(data_ft, idx_other)
    adjustment = bfft(adjustment_ft)
    adjustment .*= prod(step, freq)

    intensity = base .+ real.(adjustment)
    grid = CartesianGrid(length.(freq), ntuple(zero, length(freq)), 1 ./ step.(freq))
    georef((intensity = vec(intensity),), grid)
end

function fft_only(points::PointSet, region; nfreq, fmax)
    nufft_anydomain(
        region, nfreq, fmax, points, ones(ComplexF64, length(points)), -1, 1e-14) # TODO: make work for grid data
end
