# should be a struct implementing getargument and getestimate
# if the estimate is anisotropic, then the argument should be a tuple of length D
# if the estimate is isotropic, then the argument can either be a single `AbstractVector` or a tuple of length 1
# it is a good idea to use checkinputs when constructing these
abstract type AbstractEstimate{D, P, Q, N} end
abstract type AnisotropicEstimate{D, P, Q} <: AbstractEstimate{D, P, Q, D} end
abstract type IsotropicEstimate{D, P, Q} <: AbstractEstimate{D, P, Q, 1} end

# expects that the following two functions are implemented for any subtype
# default assumptions are that these fields exists with the names estimationinformation and processinformation
getestimationinformation(est::AbstractEstimate) = est.estimationinformation
getprocessinformation(est::AbstractEstimate) = est.processinformation
function getargument(est::AbstractEstimate)
    throw(ArgumentError("no getargument method defined for $(typeof(est))"))
end
function getestimate(est::AbstractEstimate)
    throw(ArgumentError("no getestimate method defined for $(typeof(est))"))
end
getextrainformation(est::AbstractEstimate) = () # if you need additional information, override this method
getestimatename(est::AbstractEstimate) = nameof(typeof(est)) # if you want a different name, override this method

Meshes.embeddim(::AbstractEstimate{D}) where {D} = D
Base.size(::AbstractEstimate{D, P, Q}) where {D, P, Q} = (P, Q)
function processnames(estimate::AbstractEstimate)
    processinfo = getprocessinformation(estimate)
    (processinfo.process_indices_1, processinfo.process_indices_2)
end

struct ProcessInformation{D, I1, I2, M, A}
    process_indices_1::I1
    process_indices_2::I2
    mean_product::M
    atoms::A
    function ProcessInformation(
            process_indices_1,
            process_indices_2,
            mean_product,
            atoms,
            ::Val{D}
    ) where {D}
        new{D, typeof(process_indices_1), typeof(process_indices_2),
            typeof(mean_product), typeof(atoms)}(
            process_indices_1,
            process_indices_2,
            mean_product,
            atoms
        )
    end
end
Base.ndims(::ProcessInformation{D}) where {D} = D

function getprocessinformationindex(est::AbstractEstimate, p, q)
    processinformation = getprocessinformation(est)
    return ProcessInformation(
        processinformation.process_indices_1[p],
        processinformation.process_indices_2[q],
        processinformation.mean_product[p, q],
        processinformation.atoms[p, q],
        Val{ndims(processinformation)}()
    )
end

function checkprocessinformation(
        processinformation::ProcessInformation, P, Q)
    if length(processinformation.process_indices_1) != P
        throw(ArgumentError("processinformation.process_indices_1 should have length $P"))
    end
    if length(processinformation.process_indices_2) != Q
        throw(ArgumentError("processinformation.process_indices_2 should have length $Q"))
    end
    if size(processinformation.mean_product, 1) != P ||
       size(processinformation.mean_product, 2) != Q
        throw(ArgumentError(
            "processinformation.mean_product should have size ($P, $Q), but has size $(size(processinformation.mean_product))."
        ))
    end
    if size(processinformation.atoms, 1) != P || size(processinformation.atoms, 2) != Q
        throw(ArgumentError(
            "processinformation.atoms should have size ($P, $Q), but has size $(size(processinformation.atoms))."
        ))
    end
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
        estimate::AbstractEstimate{D, P, Q, N},
        p,
        q,
        i::Vararg{Any, N}
) where {D, P, Q, N}
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
        estimate::AbstractEstimate{D, P, Q, N}, i::Vararg{Any, N}) where {D, P, Q, N}
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
function Base.getindex(
        estimate::AbstractEstimate{D, P, Q, N},
        p,
        q,
        i::Vararg{Any, N}
) where {D, Q, P, N}
    checkbounds(estimate, p, q, i...)
    return constructorof(typeof(estimate))(
        getargumentindex(estimate, i...),
        getestimateindex(estimate, p, q, i...),
        getprocessinformationindex(estimate, p, q),
        getestimationinformation(estimate),
        getextrainformation(estimate)...
    )
end

function Base.getindex(estimate::AbstractEstimate, p, q)
    checkbounds(estimate, p, q)
    return constructorof(typeof(estimate))(
        getargument(estimate),
        getestimateindex(estimate, p, q),
        getprocessinformationindex(estimate, p, q),
        getestimationinformation(estimate),
        getextrainformation(estimate)...
    )
end

## get argument index
function getargumentindex(
        estimate::AbstractEstimate{D, P, Q, N}, i::Vararg{Any, N}) where {D, P, Q, N}
    getindex.(getargument(estimate), tuple(i...))
end

## get estimate index
function getestimateindex(
        estimate::AbstractEstimate{D, P, Q, N},
        p,
        q,
        i::Vararg{Any, N}
) where {D, P, Q, N}
    _getestimateindex(getestimate(estimate), p, q, i...)
end

function getestimateindex(estimate::AbstractEstimate{D, P, Q, N}, p, q) where {D, P, Q, N}
    _getestimateindex(getestimate(estimate), p, q)
end

function _getestimateindex(
        estimate::AbstractArray{<:SMatrix, N},
        p,
        q,
        i::Vararg{Any, N}
) where {N}
    _getestimateindex(estimate, p, q)[i...]
end
function _getestimateindex(estimate::AbstractArray{<:SMatrix, N}, p, q) where {N}
    getindex.(estimate, p, q)
end

function _getestimateindex(
        estimate::AbstractArray{<:Number, M},
        p,
        q,
        i::Vararg{Any, N}
) where {M, N}
    estimate[p, q, i...]
end
function _getestimateindex(estimate::AbstractArray{<:Number, N}, p, q) where {N}
    collect(selectdim(selectdim(estimate, 1, p), 1, q)) # once first is selected, second becomes first
end

function _getestimateindex(
        estimate::AbstractArray{<:Number, N},
        p,
        q,
        i::Vararg{Any, N}
) where {N} # one process case, p and q have been checked by this point but should be 1
    estimate[i...]
end

function _getestimateindex(
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}},
        p,
        q,
        i::Vararg{Any, N}
) where {N}
    return Dict(key => val[i...] for (key, val) in _getestimateindex(estimate, p, q))
end

function _getestimateindex(
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}},
        p,
        q
) where {N}
    return Dict(
        (pᵢ, qᵢ) => estimate[(pᵢ, qᵢ)] for (pᵢ, qᵢ) in Iterators.product(p, q) if p ≤ q
    ) # assumes that Dict is only used when we have symmetric cases
end

function _getestimateindex(
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}},
        p::Int,
        q::Int
) where {N}
    return estimate[(min(p, q), max(p, q))] # when we only have one process of each we give back the values
end

## constructor input checking
function checkinputs(
        argument::NTuple{D, F},
        estimate::AbstractArray{SMatrix{P, Q, T, L}, N},
        processinformation::ProcessInformation
) where {D, F, P, Q, T, L, N}
    @assert N==D "frequencies and estimate should be same dimension when estimate is an array of matrices."
    @assert length.(argument) == size(estimate)
    checkprocessinformation(processinformation, P, Q)
    return P, Q
end

function checkinputs(
        argument::AbstractVector, estimate, processinformation::ProcessInformation)
    checkinputs((argument,), estimate, processinformation)
end

function checkinputs(argument::NTuple{D, F}, estimate::AbstractArray{<:Number, D},
        processinformation::ProcessInformation) where {D, F}
    @assert length.(argument)==size(estimate) "argument should have same length as each dimension of estimate, but got $(length.(argument)) and $(size(estimate))"
    checkprocessinformation(processinformation, 1, 1)
    return 1, 1
end

function checkinputs(
        argument::NTuple{D, F}, estimate::AbstractArray{T, N},
        processinformation::ProcessInformation) where {D, F, T, N}
    @assert N==D + 2 "argument should either have same dimension as estimate (if one process being analysed) or the dim(estimate)==dim(freq)+2."
    @assert length.(argument) == size(estimate)[3:end]
    checkprocessinformation(processinformation, size(estimate, 1), size(estimate, 2))
    return size(estimate, 1), size(estimate, 2)
end

# TODO: dictionary case doesn't currently work in all situations
function checkinputs(
        argument::NTuple{D, F},
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, D}},
        processinformation::ProcessInformation
) where {D, F}
    @assert all(x -> size(x) == length.(argument), values(estimate)) "all estimates in the dictionary should have the same size as the length.(argument)."
    P = maximum(x -> x[2], keys(estimate)) # assumes that Dict is only used when we have symmetric cases
    checkprocessinformation(processinformation, P, P)
    return P, P
end
