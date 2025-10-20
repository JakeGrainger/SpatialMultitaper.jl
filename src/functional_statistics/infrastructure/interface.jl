# should be a struct implementing getargument and getestimate
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
getestimationinformation(est::AbstractEstimate) = est.estimationinformation
getprocessinformation(est::AbstractEstimate) = est.processinformation

function is_same_process_sets(est::AbstractEstimate)
    is_same_process_sets(getprocessinformation(est))
end

function getargument(est::AbstractEstimate)
    throw(ArgumentError("no getargument method defined for $(typeof(est))"))
end
function getestimate(est::AbstractEstimate)
    throw(ArgumentError("no getestimate method defined for $(typeof(est))"))
end
getextrainformation(::AbstractEstimate) = () # if you need additional information, override this method
getbaseestimatename(est) = getbaseestimatename(typeof(est))
function getbaseestimatename(T::Type{<:AbstractEstimate}) # please override this if you want a different name, you should call `getestimatename` when using this as some types may overload that in certain cases
    typename = string(nameof(T))
    # Convert CamelCase to lowercase words with spaces
    lowercase(join(split(typename, r"(?=[A-Z])"), " "))
end
function getshortbaseestimatename(T::Type{<:AbstractEstimate})
    getbaseestimatename(T::Type{<:AbstractEstimate})
end
