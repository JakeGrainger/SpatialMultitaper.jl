"""
    rotational_spectra(args...; kwargs...)

Compute rotationally averaged spectral estimates from spatial data.

This function computes a spectral estimate and then averages it over circular annuli in
wavenumber space, providing an isotropic (rotation-invariant) summary of the wavenumber
domain variability. The result is real-valued as it is the rotationally averaged power.

# Arguments
For arguments, see [`spectra`](@ref) documentation.

# Keywords
- `radii::Union{Nothing, AbstractVector} = nothing`: Radial wavenumber values for
    averaging. If `nothing`, default radii are computed based on the wavenumber grid.
- Additional keywords: All keywords from [`spectra`](@ref) are supported.

# Returns
- `RotationalSpectra`: A rotationally averaged spectral estimate with radial wavenumber
    coordinates and real-valued power estimates.

# Examples
```julia
# Basic rotational spectrum
rot_spec = rotational_spectra(data; kmax = 0.5, nw = 3)

# With custom radii
custom_radii = 0.1:0.1:0.5
rot_spec = rotational_spectra(data, region; radii = custom_radii, kmax = 0.5)
```

See also: [`spectra`](@ref), [`rotational_estimate`](@ref)
"""
function rotational_spectra(args...; kwargs...)
    spectrum = spectra(args...; kwargs...)
    return rotational_estimate(spectrum; kwargs...)
end

"""
    rotational_partial_spectra(args...; kwargs...)

Compute rotationally averaged partial spectral estimates from spatial data.

This function computes partial spectral estimates and then averages them over circular
annuli in wavenumber space. The result is real-valued as it represents the rotationally
averaged partial power.

# Arguments
For arguments, see [`partial_spectra`](@ref) documentation.

# Keywords
- `radii::Union{Nothing, AbstractVector} = nothing`: Radial wavenumber values for
    averaging. If `nothing`, default radii are computed based on the wavenumber grid.
- Additional keywords: All keywords from [`partial_spectra`](@ref) are supported.

# Returns
- `RotationalEstimate{PartialTrait}`: A rotationally averaged partial spectral estimate
    with radial wavenumber coordinates and real-valued partial power estimates.

# Examples
```julia
# Basic rotational partial spectrum
rot_partial = rotational_partial_spectra(data; kmax = 0.5, nw = 3)

# With custom radii
custom_radii = 0.1:0.1:0.5
rot_partial = rotational_partial_spectra(data, region; radii = custom_radii, kmax = 0.5)
```

See also: [`partial_spectra`](@ref), [`rotational_estimate`](@ref)
"""
function rotational_partial_spectra(args...; kwargs...)
    spectrum = partial_spectra(args...; kwargs...)
    return rotational_estimate(spectrum; kwargs...)
end

"""
    rotational_coherence(args...; kwargs...)

Compute rotationally averaged coherence estimates from spatial data.

This function computes coherence estimates and then averages them over circular annuli
in wavenumber space, providing an isotropic summary of the coherence. The result is
real-valued as it is the rotationally averaged.

# Arguments
For arguments, see [`coherence`](@ref) documentation.

# Keywords
- `radii::Union{Nothing, AbstractVector} = nothing`: Radial wavenumber values for
    averaging. If `nothing`, default radii are computed based on the wavenumber grid.
- Additional keywords: All keywords from [`coherence`](@ref) are supported.

# Returns
- `RotationalEstimate`: A rotationally averaged coherence estimate with radial
    wavenumber coordinates and real-valued coherence estimates.

# Examples
```julia
# Basic rotational coherence
rot_coh = rotational_coherence(data; kmax = 0.5, nw = 3)

# With custom radii
custom_radii = 0.1:0.1:0.5
rot_coh = rotational_coherence(data, region; radii = custom_radii, kmax = 0.5)
```

See also: [`coherence`](@ref), [`rotational_estimate`](@ref)
"""
function rotational_coherence(args...; kwargs...)
    coh = coherence(args...; kwargs...)
    return rotational_estimate(coh; kwargs...)
end

"""
    rotational_partial_coherence(args...; kwargs...)

Compute rotationally averaged partial coherence estimates from spatial data.

This function computes partial coherence estimates and then averages them over circular
annuli in wavenumber space, providing an isotropic summary of the coherence between
processes after removing the influence of all other processes. The result is real-valued as
it is rotationally averaged.

# Arguments
For arguments, see [`partial_coherence`](@ref) documentation.

# Keywords
- `radii::Union{Nothing, AbstractVector} = nothing`: Radial wavenumber values for
    averaging. If `nothing`, default radii are computed based on the wavenumber grid.
- Additional keywords: All keywords from [`partial_coherence`](@ref) are supported.

# Returns
- `RotationalEstimate{PartialTrait}`: A rotationally averaged partial coherence estimate
    with radial wavenumber coordinates and real-valued partial coherence estimates.

# Examples
```julia
# Basic rotational partial coherence
rot_pcoh = rotational_partial_coherence(data; kmax = 0.5, nw = 3)

# With custom radii
custom_radii = 0.1:0.1:0.5
rot_pcoh = rotational_partial_coherence(data, region; radii = custom_radii, kmax = 0.5)
```

See also: [`partial_coherence`](@ref), [`rotational_estimate`](@ref)
"""
function rotational_partial_coherence(args...; kwargs...)
    coh = partial_coherence(args...; kwargs...)
    return rotational_estimate(coh; kwargs...)
end
