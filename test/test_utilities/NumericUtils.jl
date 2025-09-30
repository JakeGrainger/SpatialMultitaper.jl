module NumericUtils

"""
    slow_dft(u, f, freq, iflag)

Reference implementation of discrete Fourier transform for testing purposes.

This is a slow but accurate implementation used to verify the correctness
of faster FFT-based implementations.

# Arguments
- `u`: Spatial locations
- `f`: Function values at those locations
- `freq`: Wavenumbers to evaluate the DFT at
- `iflag`: Direction flag (≥ 0 for forward transform, < 0 for inverse)

# Returns
Vector of DFT values at the specified wavenumbers
"""
function slow_dft(u, f, freq, iflag)
    pm = iflag ≥ 0 ? 1 : -1
    return [sum(f[i] * exp(pm * 2pi * 1im * sum(u[i] .* k)) for i in eachindex(u, f))
            for k in freq]
end

export slow_dft

end # module NumericUtils
