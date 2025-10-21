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
