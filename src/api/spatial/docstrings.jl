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
partial_c_function

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

"""
    partial_k_function(data, region; kwargs...) -> KFunction{PartialTrait}
    partial_k_function(data::SpatialData; kwargs...) -> KFunction{PartialTrait}
    partial_k_function(spectrum::NormalOrRotationalSpectra; kwargs...) -> KFunction{PartialTrait}
    partial_k_function(c::CFunction{PartialTrait}) -> KFunction{PartialTrait}

Compute partial Ripley's K function from spatial data or estimates.

The partial K function represents the direct spatial clustering or regularity patterns
after removing the linear influence of all other processes. It is computed from partial
C functions using the same transformation as regular K functions, but applied to partial
estimates that show only direct relationships between processes.

# Arguments
- `data`: Spatial data for direct partial K function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `spectrum`: Spectral estimate (partial or marginal trait)
- `c::CFunction{PartialTrait}`: Partial C function estimate

# Keywords
All keywords from [`partial_c_function`](@ref) and [`k_function`](@ref) are supported:
- `radii::AbstractVector`: Distances at which to evaluate the partial K function (required)
- Plus all spectral estimation keywords: `nk`, `kmax`, `dk`, `tapers`, `nw`, `mean_method`

# Returns
- `KFunction{PartialTrait}`: A partial K function estimate with the same structure as
    regular K functions but representing direct spatial relationships.

# Notes
- Automatically converts marginal spectra/C functions to partial when needed
- Uses the same intensity normalization and volume correction as regular K functions
- Results represent direct spatial clustering after removing indirect effects
- Cannot compute from marginal C functions (must use partial C functions or convert first)

# Examples
```julia
# Direct computation from data
radii = 0.1:0.1:2.0
partial_kf = partial_k_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing partial C function
partial_cf = partial_c_function(data; radii = radii, kmax = 0.5)
partial_kf = partial_k_function(partial_cf)

# From marginal spectrum (automatically converts to partial)
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
partial_kf = partial_k_function(marginal_spec; radii = radii)

# Compare direct vs total clustering
regular_kf = k_function(data; radii = radii, kmax = 0.5)
direct_clustering = partial_kf.value .- regular_kf.value
```

See also: [`k_function`](@ref), [`partial_c_function`](@ref), [`partial_spectra`](@ref)
"""
partial_k_function

"""
    l_function(data, region; kwargs...) -> LFunction
    l_function(data::SpatialData; kwargs...) -> LFunction
    l_function(k::KFunction; kwargs...) -> LFunction
    l_function(est::AbstractEstimate; kwargs...) -> LFunction

Compute L function from spatial data, K functions, or other estimates.

The L function is a variance-stabilizing transformation of Ripley's K function defined as:
L(r) = (K(r)/V_d)^(1/d), where V_d is the volume of a unit ball in d dimensions. This
transformation makes the L function approximately linear for Poisson processes (L(r) ≈ r),
which greatly simplifies the interpretation of clustering and regularity patterns compared
to the K function.

# Arguments
- `data`: Spatial data for direct L function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `k::KFunction`: K function estimate to transform to L function
- `est::AbstractEstimate`: Any estimate that can be converted via K function

# Keywords
When computing directly from data, all keywords from [`k_function`](@ref) are supported:
- `radii::AbstractVector`: Distances at which to evaluate the L function (required)
- Plus all keywords from [`c_function`](@ref) and [`spectra`](@ref): `wavenumber_radii`,
    `rotational_method`, `nk`, `kmax`, `dk`, `tapers`, `nw`, `mean_method`

# Returns
- `LFunction`: An L function estimate containing:
  - `radii`: Distance values where L function is evaluated
  - `value`: L function values (real-valued, can be negative)
  - `processinformation`: Information about the analyzed processes
  - `estimationinformation`: Details about the estimation procedure

# Mathematical Details
The L function transformation from K function is:

L(r) = sign(K(r)) × (|K(r)|/V_d)^(1/d)

where:
- K(r) is the input K function value
- V_d is the volume of a unit d-dimensional ball:
  - V₁ = 2, V₂ = π, V₃ = 4π/3, V_d = π^(d/2)/Γ(d/2 + 1)
- The sign is preserved to handle negative K values properly
- The power 1/d provides variance stabilization

**Key Properties:**
- For Poisson processes: L(r) = r (linear baseline)
- L(r) > r indicates clustering at distance r
- L(r) < r indicates regularity/inhibition at distance r
- Better variance properties than K function for statistical testing

# Notes
- L function provides better variance stabilization than K function
- Linear behavior L(r) = r for Poisson processes simplifies interpretation
- Use [`partial_l_function`](@ref) for partial L functions from partial estimates

# Examples
```julia
# Direct computation from data and region
radii = 0.1:0.1:2.0
lf = l_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing K function
kf = k_function(data; radii = radii, kmax = 0.5)
lf = l_function(kf)

# From spectral estimate (via K function)
spec = spectra(data; kmax = 0.5, nw = 3)
lf = l_function(spec; radii = radii)
```

See also: [`k_function`](@ref), [`partial_l_function`](@ref), [`c_function`](@ref)
"""
l_function

"""
    partial_l_function(data, region; kwargs...) -> LFunction{PartialTrait}
    partial_l_function(data::SpatialData; kwargs...) -> LFunction{PartialTrait}
    partial_l_function(est::AbstractEstimate; kwargs...) -> LFunction{PartialTrait}
    partial_l_function(k::KFunction{PartialTrait}) -> LFunction{PartialTrait}

Compute partial L function from spatial data or estimates.

The partial L function represents the direct spatial clustering or regularity patterns
after removing the linear influence of all other processes. It applies the same
variance-stabilizing transformation L(r) = (K(r)/V_d)^(1/d) to partial K functions,
providing a linearized view of direct spatial relationships with the same baseline
L(r) = r for direct Poisson-like behavior.

# Arguments
- `data`: Spatial data for direct partial L function computation
- `region::Meshes.Geometry`: Spatial region for direct computation
- `est::AbstractEstimate`: Any estimate that can be converted via partial K function
- `k::KFunction{PartialTrait}`: Partial K function estimate to transform

# Keywords
All keywords from [`partial_k_function`](@ref) and [`l_function`](@ref) are supported:
- `radii::AbstractVector`: Distances at which to evaluate the partial L function (required)
- Plus all spectral estimation keywords: `nk`, `kmax`, `dk`, `tapers`, `nw`, `mean_method`

# Returns
- `LFunction{PartialTrait}`: A partial L function estimate with the same structure as
    regular L functions but representing direct spatial relationships.

# Mathematical Details
The partial L function transformation is identical to the regular L function:

L_partial(r) = sign(K_partial(r)) × (|K_partial(r)|/V_d)^(1/d)

but applied to partial K functions K_partial(r) that show only direct relationships.

**Interpretation:**
- L_partial(r) = r: no direct clustering or regularity at distance r
- L_partial(r) > r: direct clustering after removing indirect effects
- L_partial(r) < r: direct regularity/inhibition after removing indirect effects

# Notes
- Automatically converts marginal estimates to partial when needed
- Uses the same variance-stabilization properties as regular L functions
- Cannot compute from marginal K or C functions (must use partial estimates)
- Baseline comparison is still L(r) = r regardless of partial vs marginal

# Examples
```julia
# Direct computation from data
radii = 0.1:0.1:2.0
partial_lf = partial_l_function(data, region; radii = radii, kmax = 0.5, nw = 3)

# From existing partial K function
partial_kf = partial_k_function(data; radii = radii, kmax = 0.5)
partial_lf = partial_l_function(partial_kf)

# From marginal spectrum (automatically converts to partial)
marginal_spec = spectra(data; kmax = 0.5, nw = 3)
partial_lf = partial_l_function(marginal_spec; radii = radii)

```

See also: [`l_function`](@ref), [`partial_k_function`](@ref), [`partial_spectra`](@ref)
"""
partial_l_function
