"""
    partial_spectra_uncorrected(spectrum::Spectra{MarginalTrait}) -> Spectra{PartialTrait}

Compute partial spectra without finite-sample bias correction.

This function computes partial spectra using the raw inverse relationship without the
finite-sample bias correction that accounts for the number of tapers. Results may be
biased for small numbers of tapers, so use with caution. Prefer [`partial_spectra`](@ref)
for most applications.

# Arguments
- `spectrum::Spectra{MarginalTrait}`: A marginal spectral estimate

# Returns
- `Spectra{PartialTrait}`: Uncorrected partial spectral estimates

# Notes
- Equivalent to calling `partial_spectra(x, nothing)` on the spectral matrices
- Useful for theoretical analysis or when bias correction is undesired
- Results will differ from corrected partial spectra, especially with few tapers

# Examples
```julia
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
uncorrected_partial = partial_spectra_uncorrected(marginal_spec)
corrected_partial = partial_spectra(marginal_spec)
```

See also: [`partial_spectra`](@ref)
"""
function partial_spectra_uncorrected(spectrum::Spectra{MarginalTrait})::Spectra{PartialTrait}
    # Create a modified spectrum with no taper information for uncorrected computation
    new_spectrum = Spectra{MarginalTrait}(
        get_evaluation_points(spectrum), get_estimates(spectrum),
        get_process_information(spectrum), EstimationInformation(nothing))
    return partial_spectra(new_spectrum)
end
