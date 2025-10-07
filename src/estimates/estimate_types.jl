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

function Base.:(==)(a::AbstractEstimate, b::AbstractEstimate)
    getestimate(a) == getestimate(b) &&
        getargument(a) == getargument(b) &&
        getprocessinformation(a) == getprocessinformation(b) &&
        getestimationinformation(a) == getestimationinformation(b) &&
        getextrainformation(a) == getextrainformation(b)
end

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
getestimatename(T::Type{<:MarginalAbstractEstimate}) = getbaseestimatename(T)
getestimatename(T::Type{<:PartialAbstractEstimate}) = "partial " * getbaseestimatename(T)
getestimatename(est::AbstractEstimate) = getestimatename(typeof(est))

Meshes.embeddim(::AbstractEstimate{E, D}) where {E, D} = D
Base.size(est::AbstractEstimate) = length.(processnames(est))
function processnames(estimate::AbstractEstimate)
    processinfo = getprocessinformation(estimate)
    (processinfo.process_indices_1, processinfo.process_indices_2)
end

"""
    process_trait(estimate::AbstractEstimate)

Determine the trait for an estimate based on its process information.
"""
function process_trait(estimate::AbstractEstimate)
    process_info = getprocessinformation(estimate)
    _process_trait_from_info(process_info)
end

function getprocessinformationindex(est::AbstractEstimate, p, q)
    processinformation = getprocessinformation(est)
    processtrait = index_process_trait(est, p, q)
    return _getprocessinformationindex(processinformation, processtrait, p, q)
end

function _getprocessinformationindex(
        processinformation, processtrait::SingleProcessTrait, p, q)
    return ProcessInformation{embeddim(processinformation), typeof(processtrait)}(
        processinformation.process_indices_1[p],
        processinformation.process_indices_2[q],
        processinformation.mean_product[p, q],
        processinformation.atoms[p, q]
    )
end
function _getprocessinformationindex(
        processinformation, ::MultipleVectorTrait, p::Int, q::Int)
    _getprocessinformationindex(processinformation, SingleProcessTrait(), p, q)
end
function _getprocessinformationindex(
        processinformation, ::MultipleTupleTrait, p::Int, q::Int)
    _getprocessinformationindex(processinformation, SingleProcessTrait(), p, q)
end
function _getprocessinformationindex(
        processinformation, processtrait::MultipleVectorTrait, p, q)
    return ProcessInformation{embeddim(processinformation), typeof(processtrait)}(
        processinformation.process_indices_1[p],
        processinformation.process_indices_2[q],
        processinformation.mean_product[_index_to_vec(p), _index_to_vec(q)],
        processinformation.atoms[_index_to_vec(p), _index_to_vec(q)]
    )
end

function _getprocessinformationindex(
        processinformation, processtrait::MultipleTupleTrait, p, q)
    return ProcessInformation{embeddim(processinformation), typeof(processtrait)}(
        processinformation.process_indices_1[p],
        processinformation.process_indices_2[q],
        processinformation.mean_product[_index_to_svec(p), _index_to_svec(q)],
        processinformation.atoms[_index_to_svec(p), _index_to_svec(q)]
    )
end

struct EstimationInformation{T <: Union{Nothing, Int}}
    ntapers::T
end

"""
    collect(estimate::AbstractEstimate)

Produces a Tuple containing the arguments and the estimate.
"""
function Base.collect(estimate::AbstractEstimate)
    tuple(getargumenttuple(estimate)..., getestimate(estimate))
end
function getargumenttuple(estimate::AbstractEstimate)
    getargument(estimate)::Tuple
end
function getargumenttuple(estimate::IsotropicEstimate)
    _argument2tuple(getargument(estimate))
end
_argument2tuple(x::Tuple) = x
_argument2tuple(x) = (x,)

## bounds checking
function Base.checkbounds(
        estimate::AbstractEstimate{E, D, N},
        p,
        q,
        i::Vararg{Any, N}
) where {E, D, N}
    checkprocessbounds(estimate, p, q)
    checkindexbounds(estimate, i...)
    # TODO should have option to disable
    nothing
end
function Base.checkbounds(estimate::AbstractEstimate, p, q)
    checkprocessbounds(estimate, p, q)
    nothing
end

## checking process bounds
function checkprocessbounds(estimate::AbstractEstimate, p, q)
    _checkprocessbounds(getestimate(estimate), p, q)
    nothing
end

function _checkprocessbounds(estimate::AbstractArray{<:Number, M}, p, q) where {M}
    checkbounds(Bool, axes(estimate, 1), p) || throw(
        BoundsError(
        "The first process index, $p, is out of bounds for estimate with size $(size(estimate))",
    ),
    )
    checkbounds(Bool, axes(estimate, 2), q) || throw(
        BoundsError(
        "The second process index, $q, is out of bounds for estimate with size $(size(estimate))",
    ),
    )
end

function _checkprocessbounds(
        ::AbstractArray{SMatrix{P, Q, T, L}, D}, p, q) where {P, Q, T, L, D}
    checkbounds(Bool, 1:P, p) || throw(
        BoundsError(
        "The first process index, $p, is out of bounds for estimate with size $((P,Q))",
    ),
    )
    checkbounds(Bool, 1:Q, q) || throw(
        BoundsError(
        "The second process index, $q, is out of bounds for estimate with size $((P,Q))",
    ),
    )
end

function _checkprocessbounds(
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}},
        p,
        q
) where {N}
    [_checkprocessbounds(estimate, pᵢ, qᵢ) for (pᵢ, qᵢ) in Iterators.product(p, q)] # if one is an Int this still iterates, if both are Int then method below is used
end

function _checkprocessbounds(
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}},
        p::Int,
        q::Int
) where {N}
    (min(p, q), max(p, q)) in keys(estimate) || throw(
        BoundsError(
        "Process indices ($p, $q) are out of bounds for estimate with keys $(keys(estimate))",
    ),
    )
    nothing
end

## checking index bounds
function checkindexbounds(
        estimate::AbstractEstimate{E, D, N}, i::Vararg{Any, N}) where {E, D, N}
    _checkindexbounds(getestimate(estimate), i...)
end

function _checkindexbounds(
        estimate::AbstractArray{<:Number, M},
        i::Vararg{Any, N}
) where {M, N}
    checkbounds(
        Bool,
        selectdim(
            selectdim(estimate, 1, first(axes(estimate, 1))),
            1,
            first(axes(estimate, 2))
        ),
        i...
    ) || throw(
        BoundsError(
        "Index $i is out of bounds for estimate with size $(size(estimate)[3:end])",
    ),
    )
end

function _checkindexbounds(
        estimate::AbstractArray{SMatrix{P, Q, T, L}, D},
        i::Vararg{Any, D}
) where {P, Q, T, L, D}
    checkbounds(Bool, estimate, i...) || throw(
        BoundsError("Index $i is out of bounds for estimate with size $(size(estimate))"),
    )
end

function _checkindexbounds(
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}},
        i::Vararg{Any, N}
) where {N}
    if !all(x -> checkbounds(Bool, x, i...), values(estimate))
        all(size(x) == size(first(values(estimate))), values(estimate)) || throw(
            BoundsError(
            "In an estimate, all estimates in the dictionary should have the same size, these do not.",
        ),
        )
        throw(
            BoundsError(
            "Index $i is out of bounds for estimate with size $(size(first(values(estimate))))",
        ),
        )
    end
    nothing
end

## getindex
"""
    Base.getindex(estimate::AbstractEstimate, p, q)
    Base.getindex(estimate::AbstractEstimate, p, q, i...)

Get a subset of the estimate corresponding to processes `p` and `q` and potentially spatial
indices `i`.

If `p` and `q` are integers, the result will be a single process estimate.
Otherwise the result will be the same type of estimate but with the specified processes.
If one of the indices is an integer and the other is not, the result will still have the
original number of dimensions internally, in the sense that if results are stored as `Array`
then the size would be `(1, Q, ...)` or `(P, 1, ...)` as appropriate.
If they are `Array{SMatrix}` then the size of the array is the same, but each element is an
`SMatrix{1, Q}` or `SMatrix{P, 1}` as appropriate.
"""
function Base.getindex(
        estimate::AbstractEstimate{E, D, N},
        p,
        q,
        i::Vararg{Any, N}
) where {E, D, N}
    checkbounds(estimate, p, q, i...)
    return _construct_estimate_subset(
        typeof(estimate),
        E,
        getargumentindex(estimate, i...),
        getestimateindex(estimate, p, q, i...),
        getprocessinformationindex(estimate, p, q),
        getestimationinformation(estimate),
        getextrainformation(estimate)...
    )
end

function Base.getindex(estimate::AbstractEstimate{E}, p, q) where {E}
    checkbounds(estimate, p, q)
    return _construct_estimate_subset(
        typeof(estimate),
        E,
        getargument(estimate),
        getestimateindex(estimate, p, q),
        getprocessinformationindex(estimate, p, q),
        getestimationinformation(estimate),
        getextrainformation(estimate)...
    )
end

# Fallback for types that need trait parameter explicitly
function _construct_estimate_subset(
        ::Type{T}, trait::Type{<:EstimateTrait}, args...) where {T <: AbstractEstimate}
    # Extract the base constructor and call it with the trait
    base_constructor = constructorof(T)
    return base_constructor{trait}(args...)
end

## get argument index
function getargumentindex(
        estimate::AbstractEstimate{E, D, N}, i::Vararg{Any, N}) where {E, D, N}
    _getargumentindex(getargument(estimate), i...)
end
function _getargumentindex(argument::NTuple{N}, i::Vararg{Any, N}) where {N}
    getindex.(argument, tuple(i...))
end
function _getargumentindex(argument, i)
    argument[i]
end

## get estimate index
function getestimateindex(
        estimate::AbstractEstimate{E, D, N}, p, q, i::Vararg{Any, N}) where {E, D, N}
    return _getestimateindex(process_trait(estimate), getestimate(estimate), p, q, i...)
end

function getestimateindex(
        estimate::AbstractEstimate{E, D, N}, p, q) where {E, D, N}
    return _getestimateindex(process_trait(estimate), getestimate(estimate), p, q)
end

function _getestimateindex(::SingleProcessTrait, estimate::AbstractArray{<:Number, N},
        p::Int, q::Int, i::Vararg{Any, N}) where {N}
    return estimate[i...]
end
function _getestimateindex(::SingleProcessTrait, estimate::AbstractArray{<:Number},
        p::Int, q::Int)
    return estimate
end

function _getestimateindex(::MultipleTupleTrait, estimate::AbstractArray{<:SMatrix, N},
        p, q, i::Vararg{Any, N}) where {N}
    return _getestimateindex(MultipleTupleTrait(), estimate, p, q)[i...]
end
function _getestimateindex(
        ::MultipleTupleTrait, estimate::AbstractArray{<:SMatrix, N}, p, q) where {N}
    return getindex.(estimate, Ref(_index_to_svec(p)), Ref(_index_to_svec(q)))
end
function _getestimateindex(::MultipleTupleTrait, estimate::AbstractArray{<:SMatrix, N},
        p::Int, q::Int) where {N}
    return getindex.(estimate, p, q)
end

_index_to_svec(i::Int) = SVector(i)
_index_to_svec(i::SVector) = i
_index_to_svec(i::SOneTo) = i
_index_to_svec(i::AbstractVector) = SVector(i...)

function _getestimateindex(::MultipleVectorTrait, estimate::AbstractArray{<:Number, M},
        p::Int, q::Int, i::Vararg{Any, N}) where {M, N}
    estimate[p, q, i...]
end

function _getestimateindex(::MultipleVectorTrait, estimate::AbstractArray{<:Number, M},
        p, q, i::Vararg{Any, N}) where {M, N}
    estimate[_index_to_vec(p), _index_to_vec(q), i...]
end
_index_to_vec(i::Int) = i:i
_index_to_vec(i::AbstractVector) = i

function _getestimateindex(::MultipleVectorTrait, estimate::AbstractArray{<:Number, N},
        p::Int, q::Int) where {N}
    estimate[p, q, ntuple(Returns(:), Val{N - 2}())] # once first is selected, second becomes first
end
function _getestimateindex(
        ::MultipleVectorTrait, estimate::AbstractArray{<:Number, N}, p, q) where {N}
    estimate[_index_to_vec(p), _index_to_vec(q), ntuple(Returns(:), Val{N - 2}())...] # once first is selected, second becomes first
end

## constructor input checking
function checkinputs(
        argument::AbstractVector, estimate, processinformation::ProcessInformation)
    checkinputs((argument,), estimate, processinformation)
end

function checkinputs(argument::NTuple{D}, estimate::AbstractArray{T, N},
        processinformation::ProcessInformation) where {D, T, N}
    @argcheck length(argument) <= ndims(estimate) <= length(argument) + 2
    @argcheck length.(argument) == size(estimate)[(N + 1 - D):end]
    return checkprocessinformation(processinformation, estimate)
end

function checkinputs(argument::Union{NTuple{N, <:Number}, Number}, estimate::Number,
        processinformation::ProcessInformation) where {N}
    return nothing
end
