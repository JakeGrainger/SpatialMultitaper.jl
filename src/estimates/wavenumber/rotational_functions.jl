"""
    rotational_spectra(args...; radii = nothing, kwargs...)

Computes a rotationally averaged spectrum.

See `spectra` for argument details.
"""
function rotational_spectra(args...; radii = nothing, kwargs...)
    spectrum = spectra(args...; kwargs...)
    radii = _process_rotational_radii(spectrum, radii)
    real(rotational_estimate(spectrum; radii = radii)) # real valued as isotropic
end

"""
    partial_rotational_spectra(args...; radii = nothing, kwargs...)

Computes a rotationally averaged partial spectrum.

See `partial_spectra` for argument details.
"""
function partial_rotational_spectra(args...; radii = nothing, kwargs...)
    spectrum = partial_spectra(args...; kwargs...)
    radii = _process_rotational_radii(spectrum, radii)
    real(rotational_estimate(spectrum; radii = radii)) # real valued as isotropic
end

"""
    rotational_coherence(args...; radii = nothing, kwargs...)

Computes a rotationally averaged coherence.

See `spectra` for argument details.
"""
function rotational_coherence(args...; radii = nothing, kwargs...)
    coh = coherence(args...; kwargs...)
    radii = _process_rotational_radii(coh, radii)
    real(rotational_estimate(coh; radii = radii)) # real valued as isotropic
end

"""
    rotational_partial_coherence(args...; radii = nothing, kwargs...)

Computes a rotationally averaged partial coherence.

See `spectra` for argument details.
"""
function rotational_partial_coherence(args...; radii = nothing, kwargs...)
    coh = partial_coherence(args...; kwargs...)
    radii = _process_rotational_radii(coh, radii)
    real(rotational_estimate(coh; radii = radii)) # real valued as isotropic
end

_process_rotational_radii(spectrum, ::Nothing) = default_rotational_radii(spectrum)
_process_rotational_radii(spectrum, radii) = radii
