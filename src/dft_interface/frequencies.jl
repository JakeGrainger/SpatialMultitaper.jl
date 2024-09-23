"""
    choose_freq_1d(nfreq, fmax)

Returns the frequencies for our fft interfaces, with `nfreq` frequencies, and `fmax` as the fmax frequency.
"""
choose_freq_1d(nfreq, fmax) = fftshift(fftfreq(nfreq, 2fmax))

"""
    freq_downsample_startindex(nfreq,oversample)

Returns the starting index for the oversampled frequencies to downsample.
"""
function freq_downsample_startindex(nfreq,oversample)
    if iseven(nfreq) || oversample==1
        return 1
    elseif iseven(oversample)
        return oversample÷2+1
    else
        return (oversample-1)÷2+1
    end
end

"""
    freq_downsample_index(nfreq, oversample)

Returns the indices required to downsample the oversampled frequencies with a given oversampling.
This assumes that the frequencies are fftshifted.

# Arguments
- `nfreq::Int`: The number of frequencies for the desired output.
- `oversample::Int`: The oversampling factor used.
"""
function freq_downsample_index(nfreq, oversample)
    freq_downsample_startindex(nfreq, oversample):oversample:nfreq*oversample
end