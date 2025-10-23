"""
    c_function(spectrum::Spectra; radii, wavenumber_radii, rotational_method) -> CFunction
    c_function(spectrum::RotationalSpectra; radii) -> CFunction
    c_function(data, region; kwargs...) -> CFunction
    c_function(data::SpatialData; kwargs...) -> CFunction

Compute spatial C function from spectral estimates or directly from spatial data.

The C function is the reduced covariance measure evaluated on punctured balls as a
function of distance. For a stationary process at radius `r`, it represents
C(r) = C̆({x : 0 < ||x|| ≤ r}), where C̆ is the reduced covariance measure. The C function
is computed via inverse Fourier transform of the power spectral density with appropriate
spatial weighting functions that depend on the dimension.

# Arguments
- `spectrum::Spectra`: Input power spectral density estimate
- `spectrum::RotationalSpectra`: Rotationally averaged spectral estimate
- `data`: Spatial data for direct C function computation
- `region::Meshes.Geometry`: Spatial region for direct computation

# Keywords
- `radii::AbstractVector`: Distances at which to evaluate the C function (required)
- `wavenumber_radii::Union{Nothing, AbstractVector} = nothing`: Radial wavenumbers for
    rotational averaging. If `nothing`, defaults are computed from the spectrum grid.
- `rotational_method = nothing`: Smoothing kernel for rotational averaging. If `nothing`,
    uses `NoRotational()` (no smoothing).

When computing directly from data, all keywords from [`spectra`](@ref) are also supported:
- `nk`: Number of wavenumbers in each dimension
- `kmax`: Maximum wavenumber in each dimension
- `dk`: Wavenumber spacing in each dimension
- `tapers`: Taper functions to use
- `nw = 3`: Space-bandwidth product for taper generation
- `mean_method = DefaultMean()`: Method for mean estimation

# Returns
- `CFunction`: A C function estimate containing:
  - `radii`: Distance values where C function is evaluated
  - `value`: C function values (real-valued)
  - `processinformation`: Information about the analyzed processes
  - `estimationinformation`: Details about the estimation procedure

# Mathematical Details
For a stationary process with power spectral density f(k), the C function is computed as:

C(r) = ∫ f(k) W(r,k) dk

where W(r,k) is the spatial weighting function that depends on dimension D:
- **1D**: W(r,k) = 2r sinc(2r|k|)
- **2D**: W(r,k) = (r/|k|) J₁(2πr|k|) for |k| > 0, πr² for k = 0
- **3D**: W(r,k) = (r/|k|)^(3/2) J₃/₂(2πr|k|) for |k| > 0, (4π/3)r³ for k = 0

where J_ν is the Bessel function of the first kind of order ν, and sinc(x) = sin(πx)/(πx).

# Notes
- The C function measures spatial correlation as a function of distance
- Results are always real-valued due to the spatial integration
- For anisotropic spectra, rotational averaging is performed first
- Zero-atom corrections are automatically applied when available
- Use [`partial_c_function`](@ref) for partial C functions from partial spectra

# Examples
```julia
# Compute C function from existing spectral estimate
spec = spectra(data; kmax = 0.5, nw = 3)
radii = 0.1:0.1:2.0
cf = c_function(spec; radii = radii)

# Direct computation from data and region
cf = c_function(data, region; radii = 0.1:0.1:2.0, nk = (64, 64), kmax = (0.5, 0.5))

# With custom rotational averaging parameters
wavenumber_radii = 0.05:0.05:0.5
cf = c_function(spec; radii = radii, wavenumber_radii = wavenumber_radii)

# From rotational spectrum (no additional averaging needed)
rot_spec = rotational_spectra(data; kmax = 0.5, nw = 3)
cf = c_function(rot_spec; radii = radii)

# Access C function values
correlation_at_1km = cf.value[findfirst(r -> r ≈ 1.0, cf.radii)]
```

See also: [`spectra`](@ref), [`partial_c_function`](@ref), [`rotational_spectra`](@ref)
"""
c_function

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
