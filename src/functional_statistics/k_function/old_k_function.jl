
function k_function(data::SpatialData; kwargs...)::KFunction
    return k_function!(c_function(data; kwargs...))
end

function k_function(est::AbstractEstimate; kwargs...)::KFunction
    mem = deepcopy(est)
    return k_function!(mem; kwargs...)
end
# function k_function!(c::CFunction{E, D})::KFunction{E, D} where {E, D}
#     mean_prod = get_process_information(c).mean_product
#     radii = get_evaluation_points(c)
#     value = _c_to_k_transform!(
#         radii, get_estimates(c), process_trait(c), mean_prod, Val{D}())
#     processinfo = get_process_information(c)
#     estimationinfo = get_estimation_information(c)
#     return KFunction{E}(radii, value, processinfo, estimationinfo)
# end

function k_function!(spectrum::NormalOrRotationalSpectra; kwargs...)::KFunction
    return k_function!(c_function(spectrum; kwargs...))
end

# Partial K functions

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

function partial_k_function(data, region; kwargs...)::KFunction{PartialTrait}
    return partial_k_function(spatial_data(data, region); kwargs...)
end

function partial_k_function(data::SpatialData; kwargs...)::KFunction{PartialTrait}
    return k_function!(partial_c_function(data; kwargs...))
end

function partial_k_function(est::AbstractEstimate; kwargs...)::KFunction{PartialTrait}
    mem = deepcopy(est)
    return partial_k_function!(mem; kwargs...)
end

function partial_k_function!(spectrum::NormalOrRotationalSpectra{PartialTrait};
        kwargs...)::KFunction{PartialTrait}
    return k_function!(spectrum; kwargs...)
end

function partial_k_function!(spectrum::NormalOrRotationalSpectra{MarginalTrait};
        kwargs...)::KFunction{PartialTrait}
    return k_function!(partial_spectra!(spectrum); kwargs...)
end

partial_k_function!(c::CFunction{PartialTrait}) = k_function!(c)

function partial_k_function!(::CFunction{MarginalTrait})
    throw(ArgumentError(
        "Cannot compute partial K function from marginal C function. " *
        "Compute from partial spectral estimates or use partial_c_function first."
    ))
end

# Internal transformation functions
