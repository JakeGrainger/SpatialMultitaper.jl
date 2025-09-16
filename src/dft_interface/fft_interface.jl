"""
	fft_anydomain(x::Array, grid::Grid, nfreq, fmax; kwargs...)

Apply an fft to data stored in `x` given that `x` is recorded on a grid `grid`.
The output is at frequencies determined by `nfreq` and `fmax`.
The output is fftshifted so that the zero frequency is in the middle of the output array.
`x` can have more dimensions than `grid` but the first dimensions (up to the `ndims(grid)`)
must be the same as the grid.
When `x` has more dimensions than `grid`, the fft is only applied over the dimensions
corresponding to the grid.
"""
function fft_anydomain(
        x::Array,
        grid::Grid,
        nfreq::NTuple{D, Int},
        fmax::NTuple{D, Number};
        kwargs...
) where {D}
    if embeddim(grid) > ndims(x)
        throw(ArgumentError("x must have at least $(embeddim(grid)) dimensions"))
    end
    if size(x)[1:embeddim(grid)] != size(grid)
        err = ArgumentError(
            """x and grid must have the same size in the first $(embeddim(grid)),
            but size of x is $(size(x)) and size of grid is $(size(grid))"""
        )
        throw(err)
    end
    # check that fmax is a integer multiple of nyquist
    nyquist = 1 ./ (2 .* unitless_spacing(grid))
    fmaxrel = round.(Int, fmax ./ nyquist)
    all(abs.(fmaxrel .- fmax ./ nyquist) .< 1e-10) ||
        throw(ArgumentError("fmax must be a multiple of nyquist"))

    # padding
    # getting `oversample` from `choose_freq_oversample` ensures that we pad to at least size(grid)
    # also ensures that the desired frequencies can be recovered at the end
    oversample = choose_freq_oversample.(size(grid), nfreq)
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
	choose_freq_oversample(n, nfreq; maxoversample=100)

Find the smallest integer `oversample` such that `nfreq*oversample ≥ n`.

This function is needed because if `nfreq ≥ n`, we can simply pad the data to size `nfreq`.
However, if `nfreq < n`, we need to pad to a larger size to be able to recover the desired
frequencies.
"""
function choose_freq_oversample(n::Int, nfreq::Int; maxoversample = 100)
    if nfreq < 1
        err = ArgumentError("nfreq must be a positive integer, but got nfreq=$(nfreq)")
        throw(err)
    end
    if n < 1
        err = ArgumentError("n must be a positive integer, but got n=$(n)")
        throw(err)
    end

    oversample = ceil(Int, n / nfreq)
    if oversample > maxoversample
        err = ErrorException("""
            Required oversampling is $(oversample) which is greater than maxoversample=$(maxoversample).
            If you use this output frequency choice, you will have very poor performance
            computationally, probably something has been missspecified.
        """)
        throw(err)
    end
    return oversample
end

"""
	unwrap_fft_output(J::AbstractArray, fmaxrel::NTuple{D, Int}) where {D}

Copies the result of the fft to give something which goes up to some multiple of the nyquist
frequency (whilst keeping the number of output frequencies).

This is needed because if we want frequencies higher than the nyquist frequency, we can use
periodicity to recover frequencies up to `fmax = fmaxrel * nyquist`.
This works because at this point we assume the fft input is centered at zero, so the output
has this periodicity.
The final output of the calling function is adjusted appropriately for shifts in the input.
"""
function unwrap_fft_output(J::AbstractArray, fmaxrel::NTuple{D, Int}) where {D}
    if ndims(J) < D
        err = ArgumentError(
            "J must have at least $(D) dimensions, but got ndims(J)=$(ndims(J))"
        )
        throw(err)
    end
    if !all(fmaxrel .> 0)
        err = ArgumentError("fmaxrel must be a tuple of positive integers, but got fmaxrel=$(fmaxrel)")
        throw(err)
    end
    if isempty(J)
        err = ArgumentError("J must be non-empty")
        throw(err)
    end

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
	unwrap_index(n::Int, a::Int)

The indices required to unwrap the fft output to go up to `n` with a given oversampling `a`.
"""
function unwrap_index(n::Int, a::Int)
    if a < 1
        err = ArgumentError("a must be a positive integer, but got a=$(a)")
        throw(err)
    end

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
