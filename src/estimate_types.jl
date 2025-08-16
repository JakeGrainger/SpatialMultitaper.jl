# should be a struct implementing getargument and getestimate
# it is a good idea to use checkinputs when constructing these
abstract type AbstractEstimate{D,P,N} end
abstract type AnisotropicEstimate{D,P} <: AbstractEstimate{D,P,D} end
abstract type IsotropicEstimate{D,P} <: AbstractEstimate{D,P,1} end

# override this behaviour if extra fields are needed
getextrafields(::AbstractEstimate) = ()

"""
    collect(estimate::AbstractEstimate)

Produces a Tuple containing the arguments and the estimate.
"""
Base.collect(estimate::AbstractEstimate) = tuple(getargument(estimate)..., getestimate(estimate))
Base.collect(estimate::IsotropicEstimate) = tuple(getargument(estimate), getestimate(estimate))

## bounds checking
function Base.checkbounds(
    estimate::AbstractEstimate{D,P,N},
    p,
    q,
    i::Vararg{Any,N},
) where {D,P,N}
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

function _checkprocessbounds(estimate::AbstractArray{<:Number,M}, p, q) where {M}
    checkbounds(Bool, axes(estimate, 1), p) || throw(BoundsError("The first process index, $p, is out of bounds for estimate with size $(size(estimate))"))
    checkbounds(Bool, axes(estimate, 2), q) || throw(BoundsError("The second process index, $q, is out of bounds for estimate with size $(size(estimate))"))
end

function _checkprocessbounds(::AbstractArray{SMatrix{P,P,T,L},D}, p, q) where {P,T,L,D}
    checkbounds(Bool, 1:P, p) || throw(BoundsError("The first process index, $p, is out of bounds for estimate with size $((P,P))"))
    checkbounds(Bool, 1:P, q) || throw(BoundsError("The second process index, $q, is out of bounds for estimate with size $((P,P))"))
end

function _checkprocessbounds(estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}}, p, q) where {N}
    [_checkprocessbounds(estimate, pᵢ, qᵢ) for (pᵢ, qᵢ) in Iterators.product(p, q)] # if one is an Int this still iterates, if both are Int then method below is used
end

function _checkprocessbounds(estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}}, p::Int, q::Int) where {N}
    (min(p,q), max(p,q)) in keys(estimate) || throw(BoundsError("Process indices ($p, $q) are out of bounds for estimate with keys $(keys(estimate))"))
    nothing
end

## checking index bounds
function checkindexbounds(estimate::AbstractEstimate{D,P,N}, i::Vararg{Any,N}) where {D,P,N}
    _checkindexbounds(getestimate(estimate), i...)
end

function _checkindexbounds(estimate::AbstractArray{<:Number,M}, i::Vararg{Any,N}) where {M,N}
    checkbounds(Bool, selectdim(selectdim(estimate, 1, first(axes(estimate,1))), 1, first(axes(estimate,2))), i...) ||
        throw(BoundsError("Index $i is out of bounds for estimate with size $(size(estimate)[3:end])"))
end

function _checkindexbounds(estimate::AbstractArray{SMatrix{P,P,T,L},D}, i::Vararg{Any,D}) where {P,T,L,D}
    checkbounds(Bool, first(estimate), i...) ||
        throw(BoundsError("Index $i is out of bounds for estimate with size $(size(first(estimate)))"))
end

function _checkindexbounds(estimate::Dict{<:Tuple, <:AbstractArray{<:Number,N}}, i::Vararg{Any,N}) where {N}
    if !all(x->checkbounds(Bool, x, i...), values(estimate))
        all(size(x) == size(first(values(estimate))), values(estimate)) ||
            throw(BoundsError("In an estimate, all estimates in the dictionary should have the same size, these do not."))
        throw(BoundsError("Index $i is out of bounds for estimate with size $(size(first(values(estimate))))"))
    end
    nothing
end

## getindex
function Base.getindex(
    estimate::AbstractEstimate{D,P,N},
    p,
    q,
    i::Vararg{Any,N},
) where {D,P,N}
    checkbounds(estimate, p, q, i...)
    return constructorof(typeof(estimate))(
        getargumentindex(estimate, i...),
        getestimateindex(estimate, p, q, i...),
        getextrafields(estimate)...,
    )
end

function Base.getindex(estimate::AbstractEstimate, p, q)
    checkbounds(estimate, p, q)
    return constructorof(typeof(estimate))(getargument(estimate), getestimateindex(estimate, p, q), getextrafields(estimate)...)
end

## get argument index
function getargumentindex(estimate::AbstractEstimate{D,P,N}, i::Vararg{Any,N}) where {D,P,N}
    getindex.(getargument(estimate), i...)
end

## get estimate index
function getestimateindex(estimate::AbstractEstimate{D,P,N}, p, q, i::Vararg{Any,N}) where {D,P,N}
    _getestimateindex(getestimate(estimate), p, q, i...)
end

function getestimateindex(estimate::AbstractEstimate{D,P,N}, p, q) where {D,P,N}
    _getestimateindex(getestimate(estimate), p, q)
end

function _getestimateindex(
    estimate::AbstractArray{<:SMatrix,N},
    p,
    q,
    i::Vararg{Any,N},
) where {N}
    _getestimateindex(estimate, p, q)[i...]
end
function _getestimateindex(estimate::AbstractArray{<:SMatrix,N}, p, q) where {N}
    getindex.(estimate, p, q)
end

function _getestimateindex(
    estimate::AbstractArray{<:Number,M},
    p,
    q,
    i::Vararg{Any,N},
) where {M, N}
    estimate[p, q, i...]
end
function _getestimateindex(estimate::AbstractArray{<:Number,N}, p, q) where {N}
    collect(selectdim(selectdim(estimate, 1, p), 1, q)) # once first is selected, second becomes first
end

function _getestimateindex(estimate::AbstractArray{<:Number,N}, p, q, i::Vararg{Any,N}) where {N} # one process case, p and q have been checked by this point but should be 1
    estimate[i...]
end

function _getestimateindex(estimate::Dict{<:Tuple, <:AbstractArray{<:Number,N}}, p, q, i::Vararg{Any,N}) where {N}
    return Dict(key => val[i...] for (key, val) in _getestimateindex(estimate, p, q))
end

function _getestimateindex(estimate::Dict{<:Tuple, <:AbstractArray{<:Number,N}}, p, q) where {N}
    return Dict((pᵢ, qᵢ) => estimate[(pᵢ, qᵢ)] for (pᵢ, qᵢ) in Iterators.product(p, q) if p ≤ q) # assumes that Dict is only used when we have symmetric cases
end

function _getestimateindex(estimate::Dict{<:Tuple, <:AbstractArray{<:Number,N}}, p::Int, q::Int) where {N}
    return estimate[(min(p,q), max(p,q))] # when we only have one process of each we give back the values
end

## constructor input checking
function checkinputs(
    argument::NTuple{D,F},
    estimate::AbstractArray{SMatrix{P,P,T,L},N},
) where {D,F,P,T,L,N}
    @assert N == D "frequencies and estimate should be same dimension when estimate is an array of matrices."
    @assert length.(argument) == size(estimate)
    return P
end

function checkinputs(argument::AbstractVector, estimate)
    checkinputs((argument,), estimate)
end

function checkinputs(argument::NTuple{D,F}, estimate::AbstractArray{<:Number,D}) where {D,F}
    @assert length.(argument) == size(estimate)
    return 1
end

function checkinputs(argument::NTuple{D,F}, estimate::AbstractArray{T,N}) where {D,F,T,N}
    @assert N == D + 2 "argument should either have same dimension as estimate (if one process being analysed) or the dim(estimate)==dim(freq)+2."
    @assert length.(argument) == size(estimate)[3:end]
    return size(estimate, 1)
end

function checkinputs(argument::NTuple{D,F}, estimate::Dict{<:Tuple, <:AbstractArray{<:Number,D}}) where {D,F}
    @assert all(x -> size(x) == length.(argument), values(estimate)) "all estimates in the dictionary should have the same size as the length.(argument)."
    P = maximum(x->x[2], keys(estimate)) # assumes that Dict is only used when we have symmetric cases
    return P
end
