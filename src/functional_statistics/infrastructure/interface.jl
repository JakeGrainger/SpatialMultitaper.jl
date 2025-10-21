# should be a struct implementing get_evaluation_points and get_estimates
# if the estimate is anisotropic, then the argument should be a tuple of length D
# if the estimate is isotropic, then the argument can either be a single `AbstractVector` or a tuple of length 1
# it is a good idea to use checkinputs when constructing these
abstract type EstimateTrait end
struct MarginalTrait <: EstimateTrait end
struct PartialTrait <: EstimateTrait end

abstract type AbstractEstimate{E <: EstimateTrait, D, N} end
const AnisotropicEstimate{E, D} = AbstractEstimate{E, D, D}
const IsotropicEstimate{E, D} = AbstractEstimate{E, D, 1}
const MarginalAbstractEstimate{D, N} = AbstractEstimate{MarginalTrait, D, N}
const PartialAbstractEstimate{D, N} = AbstractEstimate{PartialTrait, D, N}

is_partial(::MarginalAbstractEstimate) = false
is_partial(::PartialAbstractEstimate) = true

# expects that the following two functions are implemented for any subtype
# default assumptions are that these fields exists with the names estimationinformation and processinformation
get_estimation_information(est::AbstractEstimate) = est.estimationinformation
get_process_information(est::AbstractEstimate) = est.processinformation

"""
    computed_from(::Type{<:AbstractEstimate})

Define what type of estimate or data this estimate is directly computed from.
This establishes the dependency chain for automatic computation. For example,
`computed_from(::Type{Coherence}) = Spectra` indicates that coherence is
computed directly from spectra estimates.

# Arguments
- `::Type{<:AbstractEstimate}`: The estimate type to define dependencies for

# Returns
- `Type`: The type that this estimate is directly computed from

# Examples
```julia
computed_from(::Type{Coherence}) = Spectra
computed_from(::Type{KFunction}) = CFunction
computed_from(::Type{Spectra}) = SpatialData  # Base case
```
"""
function computed_from end

"""
    allocate_estimate_memory(::Type{T}, ::Type{S}; kwargs...) where {T, S}

Allocate memory structures needed for computing estimate of type T from type S.
This function should create all necessary arrays, buffers, and intermediate
storage required for the computation.

# Arguments
- `::Type{T}`: The target estimate type to allocate memory for
- `::Type{S}`: The source type that T will be computed from
- `kwargs...`: Parameters that determine memory sizes (e.g., nk, grid dimensions)

# Returns
- Memory structure suitable for computing T from S
"""
function allocate_estimate_memory end

"""
    extract_allocation_parameters(estimate::AbstractEstimate)

Extract parameters from an existing estimate that are needed for memory allocation
of dependent estimates. For example, when computing coherence from spectra, this
extracts parameters like `nk` (number of wavenumbers) from the spectra object
rather than requiring them as keyword arguments.

# Arguments
- `estimate::AbstractEstimate`: The source estimate to extract parameters from

# Returns
- `NamedTuple`: Parameters needed for allocating memory for dependent estimates
"""
function extract_allocation_parameters end

"""
    validate_core_parameters(::Type{T}; kwargs...)

Validate that the core parameters provided by the user are valid for estimate
type T. This should check parameter types, ranges, and basic consistency before
any computation begins.

# Arguments
- `::Type{T}`: The estimate type being computed
- `kwargs...`: User-provided parameters to validate

# Throws
- `ArgumentError` or similar if parameters are invalid
"""
function validate_core_parameters end

"""
    validate_memory_compatibility(::Type{T}, mem, arg; kwargs...)

Check that the preallocated memory structure is compatible with the requested
computation parameters. This should verify that array dimensions, types, and
other memory characteristics match what's needed for the computation.

# Arguments
- `::Type{T}`: The estimate type being computed
- `mem`: The preallocated memory structure
- `arg`: The input argument/data
- `kwargs...`: Computation parameters

# Throws
- Error if memory is incompatible with the requested computation
"""
function validate_memory_compatibility end

"""
    resolve_missing_parameters(::Type{T}; kwargs...)

Resolve parameters that weren't provided by the user but can be inferred or
have computation-dependent defaults. This handles interdependencies between
parameters (e.g., if nk is provided, compute dk from it).

# Arguments
- `::Type{T}`: The estimate type being computed
- `kwargs...`: User-provided parameters

# Returns
- `NamedTuple`: Resolved parameters with missing values filled in
"""
function resolve_missing_parameters end

"""
    apply_parameter_defaults(::Type{T}, resolved_params)

Apply final default values for any parameters that are still unresolved after
the resolution phase. This provides fallback values when parameters cannot
be inferred from the computation context.

# Arguments
- `::Type{T}`: The estimate type being computed
- `resolved_params`: Parameters from the resolution phase

# Returns
- `NamedTuple`: Complete parameter set with all defaults applied
"""
function apply_parameter_defaults end

"""
    compute_estimate!(::Type{T}, mem, source; kwargs...)

Perform the actual computation of estimate type T using the preallocated memory
and source data/estimate. This is where the core mathematical computation happens.

# Arguments
- `::Type{T}`: The estimate type being computed
- `mem`: Preallocated memory structure for the computation
- `source`: Source data or estimate to compute from
- `kwargs...`: Computation parameters

# Returns
- The computed estimate of type T
"""
function compute_estimate! end

"""
    get_evaluation_points(est::AbstractEstimate)

Get the set of points at which the functional statistic was evaluated.

For anisotropic estimates, returns a tuple of length D containing the evaluation points
in each dimension. For isotropic estimates, returns either a single `AbstractVector`
or a tuple of length 1.

# Arguments
- `est::AbstractEstimate`: The estimate object

# Returns
- For anisotropic estimates: Tuple of length D with evaluation points per dimension
- For isotropic estimates: AbstractVector or tuple of length 1 with evaluation points
"""
function get_evaluation_points(est::AbstractEstimate)
    throw(ArgumentError("no get_evaluation_points method defined for $(typeof(est))"))
end

"""
    get_estimates(est::AbstractEstimate)

Get the estimated values of the functional statistic.

# Arguments
- `est::AbstractEstimate`: The estimate object

# Returns
- The estimated values as computed by the statistical method
"""
function get_estimates(est::AbstractEstimate)
    throw(ArgumentError("no get_estimates method defined for $(typeof(est))"))
end

"""
    get_extra_information(est::AbstractEstimate) -> Tuple

Get additional information from an estimate object. Default implementation returns
an empty tuple. Override this method if your estimate type needs to provide
additional information beyond the standard interface.

# Arguments
- `est::AbstractEstimate`: The estimate object

# Returns
- `Tuple`: Additional information (empty by default)
"""
get_extra_information(::AbstractEstimate) = ()
