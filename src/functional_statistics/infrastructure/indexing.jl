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
    _checkprocessbounds(getestimates(estimate), p, q)
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
    _checkindexbounds(getestimates(estimate), i...)
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
        getevaluationpointsindex(estimate, i...),
        getestimatesindex(estimate, p, q, i...),
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
        getevaluationpoints(estimate),
        getestimatesindex(estimate, p, q),
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
function getevaluationpointsindex(
        estimate::AbstractEstimate{E, D, N}, i::Vararg{Any, N}) where {E, D, N}
    _getargumentindex(getevaluationpoints(estimate), i...)
end
function _getargumentindex(argument::NTuple{N}, i::Vararg{Any, N}) where {N}
    getindex.(argument, tuple(i...))
end
function _getargumentindex(argument, i)
    argument[i]
end

## get estimate index
function getestimatesindex(
        estimate::AbstractEstimate{E, D, N}, p, q, i::Vararg{Any, N}) where {E, D, N}
    return _getestimateindex(process_trait(estimate), getestimates(estimate), p, q, i...)
end

function getestimatesindex(
        estimate::AbstractEstimate{E, D, N}, p, q) where {E, D, N}
    return _getestimateindex(process_trait(estimate), getestimates(estimate), p, q)
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
