"""
	fft_anydomain(x::Array, grid::Grid, nfreq, fmax; kwargs...)

Apply an fft to data stored in `x` given that `x` is recorded on a grid `grid`.
The output is at frequencies determined by `nfreq` and `fmax`.
Note that this will have been fftshifted so that the zero frequency is in the middle of the output array.
`x` can have more dimensions than `grid` but the first dimensions must be the same as the grid (up to the size of the grid).
"""
function fft_anydomain(
        x::Array,
        grid::Grid,
        nfreq::NTuple{D, Int},
        fmax::NTuple{D, Number};
        kwargs...
) where {D}
    embeddim(grid) ≤ ndims(x) ||
        throw(ArgumentError("x must have at least $(embeddim(grid)) dimensions"))
    size(x)[1:embeddim(grid)] == size(grid) || throw(
        ArgumentError(
        "x and grid must have the same size in overlapping dimensions (the first $(embeddim(grid))), but size of x is $(size(x)) and size of grid is $(size(grid))",
    ),
    ) # need to improve this message

    # check that fmax is a integer multiple of nyquist
    nyquist = 1 ./ (2 .* unitless_spacing(grid))
    fmaxrel = round.(Int, fmax ./ nyquist)
    all(abs.(fmaxrel .- fmax ./ nyquist) .< 1e-10) ||
        throw(ArgumentError("fmax must be a multiple of nyquist"))

    # padding
    oversample = choose_freq_res.(size(grid), nfreq) # ensures that we pad to at least size(grid), but then also so that the required frequencies are present as a subarray
    padded_x_size = if ndims(x) === embeddim(grid)
        nfreq .* oversample
    else
        (nfreq .* oversample..., size(x)[(embeddim(grid) + 1):ndims(x)]...)
    end
    padded_x = padto(x, padded_x_size)

    # compute fft over the dimensions corresponding to the grid (spatial variation)
    J = fftshift(fft(padded_x, 1:embeddim(grid); kwargs...), 1:embeddim(grid))

    # downsample to desired output frequencies
    down_ind = if ndims(x) === embeddim(grid)
        CartesianIndices(freq_downsample_index.(nfreq, oversample))
    else
        CartesianIndices((
            freq_downsample_index.(nfreq, oversample)...,
            axes(x)[(embeddim(grid) + 1):ndims(x)]...
        ))
    end
    J_downsample = J[down_ind] # this does not happen inside unwrap_fft_output as this first copies, and then downsamples, which would result in very large intermediate arrays

    # unwrap fft output (if higher than the nyquist frequency was requested this copies out to that frequency, and then downsamples to the correct number of frequencies)
    J_unwrapped = unwrap_fft_output(J_downsample, fmaxrel)

    # rescale to account for shift in grid
    # temporary reshape so that we only have one residual dimension
    J_unwrapped_reshaped = if ndims(J_unwrapped) === embeddim(grid)
        reshape(J_unwrapped, (size(J_unwrapped)..., 1))
    else
        reshape(
            J_unwrapped,
            (
                size(J_unwrapped)[1:embeddim(grid)]...,
                prod(size(J_unwrapped)[(embeddim(grid) + 1):end])
            )
        )
    end

    freq = Iterators.ProductIterator(choose_freq_1d.(nfreq, fmax))
    shift = unitless_minimum(grid) .+ unitless_spacing(grid) ./ 2
    for i in axes(J_unwrapped_reshaped, ndims(J_unwrapped_reshaped))
        J_unwrapped_slice = selectdim(J_unwrapped_reshaped, ndims(J_unwrapped_reshaped), i)
        for (j, k) in zip(eachindex(J_unwrapped_slice), freq)
            J_unwrapped_slice[j] *= exp(-2pi * 1im * sum(k .* shift))
        end
    end

    return J_unwrapped
end

"""
	choose_freq_res(n, nfreq; maxitr=100)

Find the smallest integer `oversample` such that `nfreq*oversample ≥ n`.
"""
function choose_freq_res(n, nfreq; maxitr = 100)
    oversample = 1
    while oversample < maxitr
        if nfreq * oversample ≥ n
            return oversample
        end
        oversample += 1
    end
    throw(ExceptionError("Could not find oversampling within `maxitr=$maxitr` iterations."))
end

"""
	unwrap_fft_output(J::AbstractArray, fmaxrel::NTuple{D, Int}) where {D}

Copies the result of the fft to give something which goes up to some multiple of the nyquist frequency (whilst keeping the number of output frequencies).
"""
function unwrap_fft_output(J::AbstractArray, fmaxrel::NTuple{D, Int}) where {D}
    if all(fmaxrel .== 1)
        return J
    end

    ind = if ndims(J) === D
        CartesianIndices(unwrap_index.(size(J), fmaxrel))
    else
        CartesianIndices((
            unwrap_index.(size(J)[1:D], fmaxrel)..., axes(J)[(D + 1):ndims(J)]...))
    end

    repeat_times = if ndims(J) === D # repeats fmaxrel times in transformed dimensions is fmaxrel is odd, or fmaxrel+1 times if fmaxrel is even (does not repeat in the extra dimensions)
        fmaxrel .+ iseven.(fmaxrel)
    else
        (fmaxrel .+ iseven.(fmaxrel)..., ntuple(d -> 1, ndims(J) - D)...)
    end

    repeat(J, outer = repeat_times)[ind]
end

"""
	unwrap_index(n, a)

The indices required to unwrap the fft output to go up to `n` with a given oversampling `a`.
"""
function unwrap_index(n, a)
    if iseven(n) || a == 1
        if iseven(a)
            return n ÷ 2 + 1 .+ (0:a:((n - 1) * a)) # a even n even
        else
            return 1:a:(n * a) # a odd n even or a == 1
        end
    else
        if iseven(a)
            return (a + n - 1) ÷ 2 + 1 .+ (0:a:((n - 1) * a)) # a even n odd
        else
            return (((a - 1) ÷ 2) + 1):a:(n * a) # a odd n odd
        end
    end
end
