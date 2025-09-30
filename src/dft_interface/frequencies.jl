"""
	_choose_frequencies_1d(nk, kmax)

Returns the frequencies for our fft interfaces, with `nk` frequencies, and `kmax` as the kmax frequency.
"""
_choose_frequencies_1d(nk, kmax) = fftshift(fftfreq(nk, 2kmax))

"""
	freq_downsample_startindex(nk,oversample)

Returns the starting index for the oversampled frequencies to downsample.
"""
function freq_downsample_startindex(nk, oversample)
    if iseven(nk) || oversample == 1
        return 1
    elseif iseven(oversample)
        return oversample ÷ 2 + 1
    else
        return (oversample - 1) ÷ 2 + 1
    end
end

"""
	freq_downsample_index(nk, oversample)

Returns the indices required to downsample the oversampled frequencies with a given oversampling.
This assumes that the frequencies are fftshifted.

# Arguments
- `nk::Int`: The number of frequencies for the desired output.
- `oversample::Int`: The oversampling factor used.
"""
function freq_downsample_index(nk, oversample)
    freq_downsample_startindex(nk, oversample):oversample:(nk * oversample)
end
