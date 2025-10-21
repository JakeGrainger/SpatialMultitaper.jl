"""
    k_function(data, region; kwargs...) -> KFunction
    k_function(data::SpatialData; kwargs...) -> KFunction
    k_function(c::CFunction; kwargs...) -> KFunction
    k_function(spectrum::NormalOrRotationalSpectra; kwargs...) -> KFunction

Compute Ripley's K function from spatial data, C functions, or spectral estimates.

Ripley's K function K(r) measures the expected number of points within distance r of a
typical point, normalized by the intensity. It is fundamental in spatial statistics for
detecting clustering or regularity in point patterns. The K function is derived from the
C function using the relationship K(r) = C(r)/λ² + V_d r^d, where λ is the process
intensity and V_d is the volume of a unit ball in d dimensions.

# Arguments
- `data`: Spatial data for direct K function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `c::CFunction`: C function estimate to transform to K function
- `spectrum::NormalOrRotationalSpectra`: Spectral estimate (converted via C function)

# Keywords
When computing directly from data, all keywords from [`c_function`](@ref) are supported:
- `radii::AbstractVector`: Distances at which to evaluate the K function (required)
- `wavenumber_radii`: Radial wavenumbers for rotational averaging
- `rotational_method`: Smoothing kernel for rotational averaging
- Plus all keywords from [`spectra`](@ref): `nk`, `kmax`, `dk`, `tapers`, `nw`, `mean_method`

# Returns
- `KFunction`: A K function estimate containing:
  - `radii`: Distance values where K function is evaluated
  - `value`: K function values (always real and positive)
  - `processinformation`: Information about the analyzed processes
  - `estimationinformation`: Details about the estimation procedure

# Mathematical Details
For a stationary point process with intensity λ, Ripley's K function is defined as:

K_{ij}(r) = λ⁻¹ E[number of additional `i` points within distance `r` of a typical `j` point]

The transformation from C function to K function is:
K(r) = C(r)/λ² + V_d r^d

where:
- C(r) is the corresponding C function value
- λ is the process intensity (mean product of the processes)
- V_d is the volume of a unit ball in d dimensions:
  - V₁ = 2 (length of unit interval)
  - V₂ = π (area of unit circle)
  - V₃ = 4π/3 (volume of unit sphere)
  - V_d = π^(d/2) / Γ(d/2 + 1) (general formula)

# Notes
- K function values are always real and non-negative
- For Poisson processes, K(r) = V_d r^d (no clustering or regularity)
- Values above V_d r^d indicate clustering; below indicate regularity
- The function automatically handles intensity normalization from process information
- Use [`partial_k_function`](@ref) for partial K functions from partial estimates

# Examples
```julia
# Direct computation from data and region
radii = 0.1:0.1:2.0
kf = k_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing C function
cf = c_function(data; radii = radii, kmax = 0.5)
kf = k_function(cf)

# From spectral estimate (via C function)
spec = spectra(data; kmax = 0.5, nw = 3)
kf = k_function(spec; radii = radii)

# Compare to theoretical Poisson process
theoretical_poisson = π .* radii.^2  # For 2D case
clustering_measure = kf.value .- theoretical_poisson
```

See also: [`c_function`](@ref), [`partial_k_function`](@ref), [`spectra`](@ref)
"""
k_function

function k_function(data, region; kwargs...)::KFunction
    return k_function(spatial_data(data, region); kwargs...)
end

function k_function(arg; kwargs...)
    compute(KFunction, arg; kwargs...)
end
