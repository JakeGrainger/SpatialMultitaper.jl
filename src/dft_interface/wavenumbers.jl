"""
	_choose_wavenumbers_1d(nk, kmax)

Returns the wavenumbers for our fft interfaces, with `nk` wavenumbers, and `kmax` as the kmax wavenumber.
"""
_choose_wavenumbers_1d(nk, kmax) = fftshift(fftfreq(nk, 2kmax))

"""
	wavenumber_downsample_startindex(nk,oversample)

Returns the starting index for the oversampled wavenumbers to downsample.
"""
function wavenumber_downsample_startindex(nk, oversample)
    if iseven(nk) || oversample == 1
        return 1
    elseif iseven(oversample)
        return oversample ÷ 2 + 1
    else
        return (oversample - 1) ÷ 2 + 1
    end
end

"""
	wavenumber_downsample_index(nk, oversample)

Returns the indices required to downsample the oversampled wavenumbers with a given oversampling.
This assumes that the wavenumbers are fftshifted.

# Arguments
- `nk::Int`: The number of wavenumbers for the desired output.
- `oversample::Int`: The oversampling factor used.
"""
function wavenumber_downsample_index(nk, oversample)
    wavenumber_downsample_startindex(nk, oversample):oversample:(nk * oversample)
end
