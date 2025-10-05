"""
	fft_anydomain(x::Array, grid::Grid, nk, kmax; kwargs...)

Apply an fft to data stored in `x` given that `x` is recorded on a grid `grid`.
The output is at wavenumbers determined by `nk` and `kmax`.
The output is fftshifted so that the zero wavenumber is in the middle of the output array.
`x` can have more dimensions than `grid` but the first dimensions (up to the `ndims(grid)`)
must be the same as the grid.
When `x` has more dimensions than `grid`, the fft is only applied over the dimensions
corresponding to the grid.
"""
function fft_anydomain(x::Array, grid::Grid, nk::NTuple{D, Int},
        kmax::NTuple{D, Number}; kwargs...) where {D}
    mem = preallocate_fft_anydomain(x, grid, nk, kmax)
    fft_anydomain!(mem, x, grid, nk, kmax; kwargs...)
end

function preallocate_fft_anydomain(
        x::Array, grid::Grid, nk::NTuple{D, Int}, kmax::NTuple{D, Number}) where {D}
    # check that kmax is a integer multiple of nyquist
    nyquist = 1 ./ (2 .* unitless_spacing(grid))
    kmaxrel = round.(Int, kmax ./ nyquist)
    all(abs.(kmaxrel .- kmax ./ nyquist) .< 1e-10) ||
        throw(ArgumentError("kmax must be a multiple of nyquist"))

    padded_x_size, down_ind = _padding_and_downsampling_setup(x, grid, nk)

    # allocations
    padded_x = preallocate_padto_complex(x, padded_x_size)
    J = preallocate_padto_complex(x, padded_x_size)
    J_downsample = zeros(eltype(J), size(down_ind))
    J_unwrapped = preallocated_unwrap_fft_output(J_downsample, kmaxrel)

    return (padded_x = padded_x, J = J, J_downsample = J_downsample,
        J_unwrapped = J_unwrapped, down_ind = down_ind, kmaxrel = kmaxrel)
end

function fft_anydomain!(mem, x::Array, grid::Grid, nk::NTuple{D, Int},
        kmax::NTuple{D, Number}; kwargs...) where {D}
    @argcheck embeddim(grid) ≤ ndims(x)

    @argcheck size(x)[1:embeddim(grid)] == size(grid)

    # allocations
    kmaxrel = mem.kmaxrel
    down_ind = mem.down_ind
    padded_x = mem.padded_x
    J = mem.J
    J_downsample = mem.J_downsample
    J_unwrapped = mem.J_unwrapped

    # algorithm

    # padding
    padded_x = padto!(padded_x, x)

    # compute fft over the dimensions corresponding to the grid (spatial variation)
    J_x = fft!(padded_x, 1:embeddim(grid); kwargs...)
    J = fftshift!(J, J_x, 1:embeddim(grid))

    # downsample to desired output wavenumbers
    J_downsample .= view(J, down_ind) # this does not happen inside unwrap_fft_output as this first copies, and then downsamples, which would result in very large intermediate arrays

    # unwrap fft output (if higher than the nyquist wavenumber was requested this copies out to that wavenumber, and then downsamples to the correct number of wavenumbers)
    J_unwrapped = unwrap_fft_output!(J_unwrapped, J_downsample, kmaxrel)

    # rescale to account for shift in grid
    # temporary reshape so that we only have one residual dimension
    J_unwrapped_reshaped = if ndims(J_unwrapped) === embeddim(grid)
        reshape(J_unwrapped, (size(J_unwrapped)..., 1))
    else
        reshape(J_unwrapped,
            (size(J_unwrapped)[1:embeddim(grid)]...,
                prod(size(J_unwrapped)[(embeddim(grid) + 1):end])))
    end

    wavenumber = Iterators.ProductIterator(_choose_wavenumbers_1d.(nk, kmax))
    shift = unitless_minimum(grid) .+ unitless_spacing(grid) ./ 2
    for i in axes(J_unwrapped_reshaped, ndims(J_unwrapped_reshaped))
        J_unwrapped_slice = selectdim(J_unwrapped_reshaped, ndims(J_unwrapped_reshaped), i)
        for (j, k) in zip(eachindex(J_unwrapped_slice), wavenumber)
            J_unwrapped_slice[j] *= exp(-2pi * 1im * sum(k .* shift))
        end
    end

    return J_unwrapped
end

function _padding_and_downsampling_setup(x, grid, nk)
    # getting `oversample` from `choose_wavenumber_oversample` ensures that we pad to at least size(grid)
    # also ensures that the desired wavenumbers can be recovered at the end
    oversample = choose_wavenumber_oversample.(size(grid), nk)
    pad_size = (nk .* oversample..., size(x)[(embeddim(grid) + 1):ndims(x)]...)
    down_ind = CartesianIndices(
        (wavenumber_downsample_index.(nk, oversample)...,
        axes(x)[(embeddim(grid) + 1):ndims(x)]...)
    )
    return pad_size, down_ind
end

"""
	choose_wavenumber_oversample(n, nk; maxoversample=100)

Find the smallest integer `oversample` such that `nk*oversample ≥ n`.

This function is needed because if `nk ≥ n`, we can simply pad the data to size `nk`.
However, if `nk < n`, we need to pad to a larger size to be able to recover the desired
wavenumbers.
"""
function choose_wavenumber_oversample(n::Int, nk::Int; maxoversample = 100)
    if nk < 1
        err = ArgumentError("nk must be a positive integer, but got nk=$(nk)")
        throw(err)
    end
    if n < 1
        err = ArgumentError("n must be a positive integer, but got n=$(n)")
        throw(err)
    end

    oversample = ceil(Int, n / nk)
    if oversample > maxoversample
        err = ErrorException("""
            Required oversampling is $(oversample) which is greater than maxoversample=$(maxoversample).
            If you use this output wavenumber choice, you will have very poor performance
            computationally, probably something has been missspecified.
        """)
        throw(err)
    end
    return oversample
end

"""
	unwrap_fft_output(J::AbstractArray, kmaxrel::NTuple{D, Int}) where {D}

Copies the result of the fft to give something which goes up to some multiple of the nyquist
wavenumber (whilst keeping the number of output wavenumbers).

This is needed because if we want wavenumbers higher than the nyquist wavenumber, we can use
periodicity to recover wavenumbers up to `kmax = kmaxrel * nyquist`.
This works because at this point we assume the fft input is centered at zero, so the output
has this periodicity.
The final output of the calling function is adjusted appropriately for shifts in the input.
"""
function unwrap_fft_output(J::AbstractArray, kmaxrel::NTuple{D, Int}) where {D}
    J_unwrapped = preallocated_unwrap_fft_output(J, kmaxrel)
    if all(kmaxrel .== 1)
        return J
    end
    return unwrap_fft_output!(J_unwrapped, J, kmaxrel)
end

function preallocated_unwrap_fft_output(J::AbstractArray, kmaxrel::NTuple{D, Int}) where {D}
    @argcheck ndims(J) ≥ D
    @argcheck all(kmaxrel .> 0)
    @argcheck !isempty(J)

    idxs = _make_unwrap_index(J, kmaxrel)
    return zeros(eltype(J), size(idxs))
end

function unwrap_fft_output!(J_unwrapped::AbstractArray, J::AbstractArray, kmaxrel)
    idxs = _make_unwrap_index(J, kmaxrel)
    J_circular = CircularArray(J)
    for (i, idx) in enumerate(idxs)
        J_unwrapped[i] = J_circular[idx]
    end
    return J_unwrapped
end

function _make_unwrap_index(J::AbstractArray{T, D}, kmaxrel::NTuple{D, Int}) where {D, T}
    return CartesianIndices(unwrap_index.(size(J), kmaxrel))
end

function _make_unwrap_index(J::AbstractArray{T, N}, kmaxrel::NTuple{D, Int}) where {D, T, N}
    @argcheck N > D

    return CartesianIndices((
        unwrap_index.(size(J)[1:D], kmaxrel)..., axes(J)[(D + 1):N]...))
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
