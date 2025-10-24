# should be a struct implementing get_evaluation_points and get_estimates
# if the estimate is anisotropic, then the argument should be a tuple of length D
# if the estimate is isotropic, then the argument can either be a single `AbstractVector` or a tuple of length 1
# it is a good idea to use checkinputs when constructing these
# you will need to define the following:
# function computed_from end
# function allocate_estimate_memory end
# function extract_relevant_memory end
# function validate_core_parameters end
# function resolve_missing_parameters end
# function validate_memory_compatibility end
# function compute_estimate! end
# function get_evaluation_points end
# function get_estimates end
#
# see also the docstrings below

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
    allocate_estimate_memory(::Type{T}, ::Type{S}, relevant_memory; kwargs...) where {T, S}

Allocate memory structures needed for computing estimate of type T from type S.
This function should create all necessary arrays, buffers, and intermediate
storage required for the computation.

# Arguments
- `::Type{T}`: The target estimate type to allocate memory for
- `::Type{S}`: The source type that T will be computed from
- `relevant_memory`: Relevant memory structure from the previous computation step
- `kwargs...`: Parameters that determine memory sizes (e.g., nk, grid dimensions)

# Returns
- output_memory: The memory structure that will hold the computed estimate of type T
- internal_memory: Any additional internal memory needed for the computation
"""
function allocate_estimate_memory end

"""
    extract_relevant_memory(::Type{T}, source)

Extract memory or information from a source (estimate object or EstimateMemory)
that is relevant for allocating memory to compute estimate type T. This function
provides a unified interface for both computation chains and direct computation
from existing estimates.

Different estimate types may need different information:
- Some need full memory structures for reuse
- Some need only metadata (e.g., array dimensions, number of processes)
- Some need specific parameters (e.g., radii, wavenumbers)
- Some don't need any information from the source

# Arguments
- `::Type{T}`: The target estimate type that will be computed
- `source`: Either an `AbstractEstimate` object or an `EstimateMemory` structure

# Returns
- Memory structure, metadata, or other information needed by `allocate_estimate_memory`
- Return type depends on what T needs from the source
- Can return `nothing` if no information is needed

# Examples
```julia
# Extract from existing estimate object
extract_relevant_memory(::Type{CFunction}, est::Spectra) = get_estimates(est)

# Extract from EstimateMemory in computation chain
extract_relevant_memory(::Type{CFunction}, mem::EstimateMemory{<:Spectra}) = mem.power

# Extract only specific parameters
extract_relevant_memory(::Type{KFunction}, c_func::CFunction) = c_func.radii

# No extraction needed
extract_relevant_memory(::Type{SomeEstimate}, ::SpatialData) = nothing
```

# Notes
- This function handles both direct computation from estimates and computation chains
- Provides a clean, unified interface instead of separate extraction mechanisms
- The returned value is passed to `allocate_estimate_memory` as the third argument
- Should define methods for both `AbstractEstimate` and `EstimateMemory` sources when needed
"""
function extract_relevant_memory end

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
    resolve_missing_parameters(::Type{T}, arg; kwargs...)

Resolve missing parameters and handle parameter constraints specific to estimate type T.
This function handles both parameter interdependencies and default values in a single
step, allowing each estimate type to implement its own resolution logic.

Note that core parameter validation has already been performed by
`validate_core_parameters`, so this function can assume that provided parameters
are valid.

For example:
- `Spectra`: Handle `dk`/`nk`/`kmax` constraint resolution where any two determine the third
- `CFunction`: Apply defaults for `radii` based on spatial domain
- Other estimates: Implement their own parameter relationships and defaults

# Arguments
- `::Type{T}`: The estimate type being computed
- `arg`: The input argument/data (can be used to derive defaults)
- `kwargs...`: User-provided parameters (already validated by `validate_core_parameters`)

# Returns
- `NamedTuple`: Complete parameter set with all missing values resolved and constraints satisfied

# Examples
```julia
# Spectra example - handles wavenumber constraints
resolve_missing_parameters(::Type{Spectra}, data; nk=100, kwargs...)
# → Uses default dk from data, derives kmax from nk and dk

# CFunction example - applies spatial defaults
resolve_missing_parameters(::Type{CFunction}, data; kwargs...)
# → Uses default radii based on spatial domain if not provided
```
"""
function resolve_missing_parameters end

"""
    validate_memory_compatibility(::Type{T}, mem, arg; kwargs...)

Check that the preallocated memory structure is compatible with the requested
computation parameters. This should verify that array dimensions, types, and
other memory characteristics match what's needed for the computation.

# Arguments
- `::Type{T}`: The estimate type being computed
- `mem`: The preallocated memory structure (which is not the `EstimateMemory` itself, but rather the internal data it holds)
- `arg`: The input argument/data
- `kwargs...`: Computation parameters

# Throws
- Error if memory is incompatible with the requested computation
"""
function validate_memory_compatibility end

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
function get_evaluation_points end

"""
    get_estimates(est::AbstractEstimate)

Get the estimated values of the functional statistic.

# Arguments
- `est::AbstractEstimate`: The estimate object

# Returns
- The estimated values as computed by the statistical method
"""
function get_estimates end

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
