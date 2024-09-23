"""
    nufft_anydomain(region, nfreq, fmax, points::PointSet, c, iflag, eps; kwargs...)

Computes the approximate dft of a function using the nufft algorithm, with points in any region, at frequencies defined by `fmax` and `nfreq`.
Currently only 1d, 2d and 3d are supported.
"""
function nufft_anydomain(region, nfreq, fmax, points::PointSet, c, iflag, eps; kwargs...)
    @assert length(nfreq) == length(fmax) == embeddim(points) "nfreq and fmax should have the same length as the number of dimensions of the points"

    points_coords = points2coords(points)
    if embeddim(points) === 1
        return nufft1d1_anydomain(box2sides(boundingbox(region))[1], nfreq[1], fmax[1], points_coords[1], c, iflag, eps; kwargs...)
    elseif embeddim(points) === 2
        return nufft2d1_anydomain(box2sides(boundingbox(region)), nfreq, fmax, points_coords[1], points_coords[2], c, iflag, eps; kwargs...)
    elseif embeddim(points) === 3
        return nufft3d1_anydomain(box2sides(boundingbox(region)), nfreq, fmax, points_coords[1], points_coords[2], points_coords[3], c, iflag, eps; kwargs...)
    else
        error("Only 1d, 2d and 3d are supported")
    end
end

## one dimensional case
"""
    nufft1d1_anydomain(interval, nfreq, fmax, xj, cj, iflag, eps; kwargs...)

Computes the approximate dft of a 1d function using the nufft algorithm, with points in any interval, at frequencies defined by `fmax` and `nfreq`.
`xj` are the points coordinates and `cj` are the function values.
Set `iflag<0` for -i in exponent.
"""
function nufft1d1_anydomain(interval, nfreq, fmax, xj, cj, iflag, eps; kwargs...)
    output_storage, input_data = nufft1d1_anydomain_precomp(interval, nfreq, fmax, xj, cj)
    return nufft1d1_anydomain!(output_storage, input_data, nfreq, fmax, iflag, eps; kwargs...)
end

"""
    nufft1d1_anydomain!(output_storage, input_data, nfreq, fmax, iflag, eps; kwargs...)

In place version of `nufft1d1_anydomain`. Use `nufft1d1_anydomain_precomp` to precompute the input data and memory required.
"""
function nufft1d1_anydomain!(output_storage, input_data, nfreq, fmax, iflag, eps; kwargs...)
    # unpack
    xj_rescaled, oversample_x, shift_x, cj = input_data.xj_rescaled, input_data.oversample_x, input_data.shift_x, input_data.cj
    oversampled_out, out, phase_correction = output_storage.oversampled_out, output_storage.out, output_storage.phase_correction

    # compute frequencies and downsampling
    freq_x = choose_freq_1d(nfreq, fmax)
    downsample_x = freq_downsample_index(nfreq, oversample_x)
    
    # compute
    nufft1d1!(xj_rescaled, cj, iflag, eps, oversampled_out; kwargs...)
    exp_scaling = iflag < 0 ? -2π*1im : 2π*1im
    phase_correction .= exp.(exp_scaling.*shift_x.*freq_x)
    
    # rescaling
    out .= oversampled_out[downsample_x,:] .* phase_correction
    return out
end

"""
    nufft1d1_anydomain_precomp(interval, nfreq, fmax, xj, cj)

Precomputes the input data and memory required for `nufft1d1_anydomain!`.
"""
function nufft1d1_anydomain_precomp(interval, nfreq, fmax, xj, cj)
    @assert length(cj) % length(xj) == 0 "length(cj) must be a multiple of length(xj)"
    # rescale data
    xj_rescaled, oversample_x, shift_x = rescale_points(xj, nfreq, fmax, interval)
    
    # preallocate storage
    n_transforms = length(cj) ÷ length(xj)
    oversampled_out = Array{complex(eltype(cj)),2}(undef, nfreq*oversample_x, n_transforms)
    out = Array{complex(eltype(cj)),2}(undef, nfreq, n_transforms)
    phase_correction = Vector{complex(eltype(cj))}(undef, nfreq)

    # format return
    input_data = (xj_rescaled = xj_rescaled, oversample_x = oversample_x, shift_x = shift_x, cj = cj)
    output_storage = (oversampled_out = oversampled_out, out = out, phase_correction = phase_correction)
    return output_storage, input_data
end

## two dimensional case
"""
    nufft2d1_anydomain(box, nfreq, fmax, xj, yj, cj, iflag, eps; kwargs...)

Computes the approximate dft of a 2d function using the nufft algorithm, with points in any box, at frequencies defined by `fmax` and `nfreq`.
`box` here should be a tuple of intervals, one for each dimension.
Similarly for `fmax` and `nfreq`.
`xj`, `yj` are the points coordinates and `cj` are the function values.
Set `iflag<0` for -i in exponent.
"""
function nufft2d1_anydomain(box, nfreq, fmax, xj, yj, cj, iflag, eps; kwargs...)
    output_storage, input_data = nufft2d1_anydomain_precomp(box, nfreq, fmax, xj, yj, cj)
    return nufft2d1_anydomain!(output_storage, input_data, nfreq, fmax, iflag, eps; kwargs...)
end

"""
    nufft2d1_anydomain!(output_storage, input_data, nfreq, fmax, iflag, eps; kwargs...)

In place version of `nufft2d1_anydomain`. Use `nufft2d1_anydomain_precomp` to precompute the input data and memory required.
"""
function nufft2d1_anydomain!(output_storage, input_data, nfreq, fmax, iflag, eps; kwargs...)
    # unpack
    xj_rescaled, oversample_x, shift_x, yj_rescaled, oversample_y, shift_y, cj = input_data.xj_rescaled, input_data.oversample_x, input_data.shift_x, input_data.yj_rescaled, input_data.oversample_y, input_data.shift_y, input_data.cj
    oversampled_out, out, phase_correction = output_storage.oversampled_out, output_storage.out, output_storage.phase_correction

    # compute frequencies and downsampling
    freq_x = choose_freq_1d(nfreq[1], fmax[1])
    downsample_x = freq_downsample_index(nfreq[1], oversample_x)
    freq_y = choose_freq_1d(nfreq[2], fmax[2])
    downsample_y = freq_downsample_index(nfreq[2], oversample_y)
    
    # compute
    nufft2d1!(xj_rescaled, yj_rescaled, cj, iflag, eps, oversampled_out; kwargs...)
    exp_scaling = iflag < 0 ? -2π*1im : 2π*1im
    phase_correction .= exp.(exp_scaling.*(shift_x.*freq_x .+ shift_y.*freq_y'))
    for i in axes(phase_correction, 2)
        phase_correction[:,i] .= exp.(exp_scaling.*(shift_x.*freq_x .+ (shift_y*freq_y[i])))
    end
    
    # rescaling
    out .= oversampled_out[downsample_x,downsample_y,:] .* phase_correction
    return out
end

"""
    nufft2d1_anydomain_precomp(box, nfreq, fmax, xj, yj, cj)

Precomputes the input data and memory required for `nufft2d1_anydomain!`.
"""
function nufft2d1_anydomain_precomp(box, nfreq, fmax, xj, yj, cj)
    @assert length(cj) % length(xj) == 0 "length(cj) must be a multiple of length(xj)"
    # rescale data
    xj_rescaled, oversample_x, shift_x = rescale_points(xj, nfreq[1], fmax[1], box[1])
    yj_rescaled, oversample_y, shift_y = rescale_points(yj, nfreq[2], fmax[2], box[2])
    
    # preallocate storage
    n_transforms = length(cj)÷length(xj)
    oversampled_out = Array{complex(eltype(cj)),3}(undef, nfreq[1]*oversample_x, nfreq[2]*oversample_y, n_transforms)
    out = Array{complex(eltype(cj)),3}(undef, nfreq[1], nfreq[2], n_transforms)
    phase_correction = Array{complex(eltype(cj)),2}(undef, nfreq[1], nfreq[2])

    # format return
    input_data = (xj_rescaled = xj_rescaled, oversample_x = oversample_x, shift_x = shift_x, yj_rescaled = yj_rescaled, oversample_y = oversample_y, shift_y = shift_y, cj = cj)
    output_storage = (oversampled_out = oversampled_out, out = out, phase_correction = phase_correction)
    return output_storage, input_data
end

"""
    nufft3d1_anydomain(cube, nfreq, fmax, xj, yj, zj, cj, iflag, eps; kwargs...)

Computes the approximate dft of a 3d function using the nufft algorithm, with points in any cube, at frequencies defined by `fmax` and `nfreq`.
`cube` here should be a tuple of intervals, one for each dimension.
Similarly for `fmax` and `nfreq`.
`xj`, `yj`, `zj` are the points coordinates and `cj` are the function values.
Set `iflag<0` for -i in exponent.
"""
function nufft3d1_anydomain(cube, nfreq, fmax, xj, yj, zj, cj, iflag, eps; kwargs...)
    output_storage, input_data = nufft3d1_anydomain_precomp(cube, nfreq, fmax, xj, yj, zj, cj)
    nufft3d1_anydomain!(output_storage, input_data, nfreq, fmax, iflag, eps; kwargs...)
    return output_storage.out
end

"""
    nufft3d1_anydomain!(output_storage, input_data, nfreq, fmax, iflag, eps; kwargs...)

In place version of `nufft3d1_anydomain`. Use `nufft3d1_anydomain_precomp` to precompute the input data and memory required.
"""
function nufft3d1_anydomain!(output_storage, input_data, nfreq, fmax, iflag, eps; kwargs...)
    # unpack
    xj_rescaled, oversample_x, shift_x, yj_rescaled, oversample_y, shift_y, zj_rescaled, oversample_z, shift_z, cj = input_data.xj_rescaled, input_data.oversample_x, input_data.shift_x, input_data.yj_rescaled, input_data.oversample_y, input_data.shift_y, input_data.zj_rescaled, input_data.oversample_z, input_data.shift_z, input_data.cj
    oversampled_out, out, phase_correction = output_storage.oversampled_out, output_storage.out, output_storage.phase_correction

    # compute frequencies and downsampling
    freq_x = choose_freq_1d(nfreq[1], fmax[1])
    downsample_x = freq_downsample_index(nfreq[1], oversample_x)
    freq_y = choose_freq_1d(nfreq[2], fmax[2])
    downsample_y = freq_downsample_index(nfreq[2], oversample_y)
    freq_z = choose_freq_1d(nfreq[3], fmax[3])
    downsample_z = freq_downsample_index(nfreq[3], oversample_z)

    # compute
    nufft3d1!(xj_rescaled, yj_rescaled, zj_rescaled, cj, iflag, eps, oversampled_out; kwargs...)
    exp_scaling = iflag < 0 ? -2π*1im : 2π*1im
    for i in axes(phase_correction, 3), j in axes(phase_correction, 2)
        phase_correction[:,j,i] .= exp.(exp_scaling.*(shift_x.*freq_x .+ (shift_y*freq_y[j] + shift_z*freq_z[i])))
    end

    # rescaling
    out .= oversampled_out[downsample_x,downsample_y,downsample_z,:] .* phase_correction
    return out
end

"""
    nufft3d1_anydomain_precomp(cube, nfreq, fmax, xj, yj, zj, cj)

Precomputes the input data and memory required for `nufft3d1_anydomain!`.
"""
function nufft3d1_anydomain_precomp(cube, nfreq, fmax, xj, yj, zj, cj)
    @assert length(cj) % length(xj) == 0 "length(cj) must be a multiple of length(xj)"
    # rescale data
    xj_rescaled, oversample_x, shift_x = rescale_points(xj, nfreq[1], fmax[1], cube[1])
    yj_rescaled, oversample_y, shift_y = rescale_points(yj, nfreq[2], fmax[2], cube[2])
    zj_rescaled, oversample_z, shift_z = rescale_points(zj, nfreq[3], fmax[3], cube[3])

    
    # preallocate storage
    n_transforms = length(cj) ÷ length(xj)
    oversampled_out = Array{complex(eltype(cj)),4}(undef, nfreq[1]*oversample_x, nfreq[2]*oversample_y, nfreq[3]*oversample_z, n_transforms)
    out = Array{complex(eltype(cj)),4}(undef, nfreq[1], nfreq[2], nfreq[3], n_transforms)
    phase_correction = Array{complex(eltype(cj)),3}(undef, nfreq[1], nfreq[2], nfreq[3])

    # format return
    input_data = (xj_rescaled = xj_rescaled, oversample_x = oversample_x, shift_x = shift_x, yj_rescaled = yj_rescaled, oversample_y = oversample_y, shift_y = shift_y, zj_rescaled = zj_rescaled, oversample_z = oversample_z, shift_z = shift_z, cj = cj)
    output_storage = (oversampled_out = oversampled_out, out = out, phase_correction = phase_correction)
    return output_storage, input_data
end


"""
    rescale_points(x, nfreq, fmax, interval; max_oversample = 100)

Rescales points to be in the interval [-π, π] and returns the oversampling factor.
Oversampling factor is the smallest integer such that the interval is rescaled to fit in [-π, π],
and the corresponding fmax frequency is a multiple of `fmax`.
Interval just needs to have `minimum` and `maximum` defined.
"""
function rescale_points(x, nfreq, fmax, interval; max_oversample = 100)
    side_length = maximum(interval) - minimum(interval)
    l = (nfreq)/(2fmax)
    oversample = 1
    for c ∈ 1:max_oversample
        oversample = c
        if side_length/(l*c) ≤ 1
            break
        end
    end
    oversample < max_oversample || error("Could not find oversampling factor, try increasing max_oversample")
    shift = minimum(interval)+side_length/2
    x_rescaled = (x .- shift)./(l*oversample) .* 2π
    return x_rescaled, oversample, shift
end
