# Public API

function c_function(data, region; kwargs...)::CFunction
    return c_function(spatial_data(data, region); kwargs...)
end

function c_function(data::SpatialData; radii = nothing,
        rotational_wavenumber_radii = nothing,
        rotational_method = nothing, kwargs...)::CFunction
    spectrum = spectra(data; kwargs...)
    wavenumber_radii_processed = process_c_rotational_radii(
        spectrum, rotational_wavenumber_radii)
    rotational_method_processed = process_c_rotational_kernel(spectrum, rotational_method)
    radii = _validate_radii(radii, data)
    return c_function(
        spectrum; radii = radii, wavenumber_radii = wavenumber_radii_processed,
        rotational_method = rotational_method_processed)
end
process_c_rotational_kernel(spectrum, ::Nothing) = default_c_rotational_kernel(spectrum)
process_c_rotational_kernel(spectrum, kernel) = kernel

process_c_rotational_radii(spectrum, ::Nothing) = default_c_rotational_radii(spectrum)
process_c_rotational_radii(spectrum, radii) = radii

default_c_rotational_kernel(spectrum) = NoRotational()
default_c_rotational_radii(spectrum) = default_rotational_radii(spectrum)

function c_function(spectrum::Spectra; radii,
        wavenumber_radii = process_c_rotational_kernel(spectrum, nothing),
        rotational_method = process_c_rotational_kernel(spectrum, nothing))
    rot_spec = rotational_estimate(
        spectrum, radii = wavenumber_radii, kernel = rotational_method)
    return _c_function(rot_spec, radii)
end

function c_function(spectrum::RotationalSpectra; radii)
    _c_function(spectrum, radii)
end

"""
    partial_c_function(data, region; kwargs...) -> CFunction{PartialTrait}
    partial_c_function(data::SpatialData; kwargs...) -> CFunction{PartialTrait}
    partial_c_function(spectrum::NormalOrRotationalSpectra{PartialTrait}; kwargs...) -> CFunction{PartialTrait}
    partial_c_function(spectrum::NormalOrRotationalSpectra{MarginalTrait}; kwargs...) -> CFunction{PartialTrait}

Compute spatial partial C function from spectral estimates or directly from spatial data.

The partial C function represents the direct spatial correlation structure after removing
the linear influence of all other processes. It is computed from partial spectral estimates
using the same inverse Fourier transform approach as the regular C function, but applied
to partial spectra rather than marginal spectra.

# Arguments
- `data`: Spatial data for direct partial C function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `spectrum`: Spectral estimate (partial or marginal trait)

# Keywords
- `radii::AbstractVector`: Distances at which to evaluate the partial C function (required)
- Additional keywords: All keywords from [`partial_spectra`](@ref) and [`c_function`](@ref)

# Returns
- `CFunction{PartialTrait}`: A partial C function estimate with the same structure as
    regular C functions but representing direct spatial relationships.

# Notes
- Automatically converts marginal spectra to partial spectra when needed
- Uses the same spatial weighting functions as regular C functions
- Results represent direct spatial correlations after removing indirect effects
- For partial trait inputs, applies C function transformation directly

# Examples
```julia
# Direct computation from data
radii = 0.1:0.1:2.0
partial_cf = partial_c_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing partial spectral estimate
partial_spec = partial_spectra(data; kmax = 0.5, nw = 3)
partial_cf = partial_c_function(partial_spec; radii = radii)

# From marginal spectrum (automatically converts to partial)
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
partial_cf = partial_c_function(marginal_spec; radii = radii)
```

See also: [`c_function`](@ref), [`partial_spectra`](@ref)
"""
function partial_c_function(data, region; kwargs...)::CFunction{PartialTrait}
    return partial_c_function(spatial_data(data, region); kwargs...)
end

function partial_c_function(
        data::SpatialData; radii = nothing, spectra_kwargs...)::CFunction{PartialTrait}
    radii = _validate_radii(radii, data)
    f_mt = partial_spectra(data; spectra_kwargs...)
    return c_function(f_mt; radii = radii)
end

function partial_c_function(
        spectrum::NormalOrRotationalSpectra{PartialTrait}; kwargs...
)::CFunction{PartialTrait}
    return c_function(spectrum; kwargs...)
end

function partial_c_function(
        spectrum::NormalOrRotationalSpectra{MarginalTrait}; kwargs...
)::CFunction{PartialTrait}
    return c_function(partial_spectra(spectrum); kwargs...)
end

# Internal implementation

function _c_function(spectrum::NormalOrRotationalSpectra{E}, radii) where {E}
    out = preallocate_spatial_output(spectrum, radii)
    store = precompute_c_weights(spectrum, radii)

    value = _sdf2C!(out, store, spectrum, radii)
    return CFunction{E}(
        radii, value, get_process_information(spectrum), get_estimation_information(spectrum))
end

# Core transformation functions
