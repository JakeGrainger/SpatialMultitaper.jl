
# """
#     create_intensities(data::NTuple{P,PointSet}, region; tapers, nk, kmax, radii, grid, smooth_width) where {P}

# Constructs the internal intensities for each of the processes partial on all the others.

# ``\\lambda_{X\\cdot Z}(u) = \\lambda_X - \\sum_{p} \\lambda_{Z_p} \\int \\psi_j(x)d x + \\sum_{p} \\sum_{x \\in Z_p} \\psi_j(u - x)``

# where ``\\tilde\\psi_j(k) = [f_{XZ}(k)f_{ZZ}(k)^{-1}]_j`` is the Fourier transform of the jth prediction kernel.
# """
# function create_intensities(
#         data::NTuple{P, PointSet},
#         region;
#         tapers,
#         nk,
#         kmax,
#         radii,
#         grid,
#         mean_method::MeanEstimationMethod = DefaultMean(),
#         smooth_width = nothing
# ) where {P}
#     spec = spectra(data, region; tapers = tapers, nk = nk, kmax = kmax)
#     intensity = mean_estimate(data, region, mean_method)
#     kernels = prediction_kernel(spec, radii = radii, smooth_width = smooth_width)
#     kernel_integral = integrate_prediction_kernel.(
#         Ref(radii), kernels.kernels, Val{embeddim(region)}())
#     kernels_interp = ntuple(
#         j -> ntuple(
#             p -> linear_interpolation(
#                 radii,
#                 getindex.(kernels.kernels[j], p),
#                 extrapolation_bc = 0.0
#             ),
#             Val{P - 1}()
#         ),
#         Val{P}()
#     )
#     ntuple(
#         j -> create_single_intensity(
#             j,
#             intensity,
#             kernel_integral[j],
#             kernels_interp[j],
#             data,
#             grid
#         ),
#         Val{P}()
#     )
# end

# function create_single_intensity(
#         idx,
#         intensities,
#         kernel_integral,
#         kernels_interp,
#         data,
#         grid
# )
#     additional_processes = Not(idx)
#     base_intensity = intensities[idx] -
#                      sum(kernel_integral .* collect(intensities)[additional_processes])
#     data_dep = collect(data)[additional_processes]
#     intensity = [max( # max of intensity and zero to avoid negative intensities
#                      base_intensity + sum(
#                          sum(kernels_interp[j](norm(centroid(grid, i) - x).val)
#                          for x in data_dep[j]) for j in eachindex(data_dep)
#                      ),
#                      zero(eltype(base_intensity))
#                  ) for i in eachindex(grid)]
#     return georef((intensity = intensity,), grid)
# end

function integrate_prediction_kernel(radii, kernel, ::Val{D}) where {D}
    A = 2 * pi^(D / 2) / gamma(D / 2)
    return A * step(radii) * sum(k * r^(D - 1) for (k, r) in zip(kernel, radii))
end

function prediction_kernel(spec::Spectra; radii, smooth_width = nothing)
    kernel_ft = prediction_kernel_ft(spec)
    kernels = _ft2kernel.(Ref(kernel_ft.freq), kernel_ft.kernels, Ref(radii), smooth_width)
    return (radii = radii, kernels = kernels)
end

function _ft2kernel(freq::NTuple{D}, kernel_ft, radii, smooth_width::Int) where {D}
    raw = _ft2kernel(freq, kernel_ft, radii, nothing)
    return movingaverage(raw[1], smooth_width)
end

function _ft2kernel(
        freq::NTuple{D}, kernel_ft, radii, smooth_width::Nothing = nothing) where {D}
    return [prod(step, freq) * real(
                sum(
                f * pcf_weight(radius, k, Val{D}())
            for
            (f, k) in zip(kernel_ft, Iterators.product(freq...))
            ),
            ) for radius in radii]
end

function movingaverage(x, width)
    y = similar(x)
    for i in eachindex(x)
        lo = max(1, i - div(width, 2))
        hi = min(length(x), i + div(width, 2))
        y[i] = zero(eltype(x))
        for idx in lo:hi
            y[i] += x[idx]
        end
        y[i] /= length(lo:hi)
    end
    return y
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
