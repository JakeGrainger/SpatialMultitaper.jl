# Bounds checking
function Base.checkbounds(
        estimate::AbstractEstimate{E, D, N}, p, q, i::Vararg{Any, N}) where {E, D, N}
    check_process_bounds(estimate, p, q)
    check_index_bounds(estimate, i...)
    # TODO: should have option to disable
    nothing
end

function Base.checkbounds(estimate::AbstractEstimate, p, q)
    check_process_bounds(estimate, p, q)
    nothing
end

# Process bounds checking
function check_process_bounds(estimate::AbstractEstimate, p, q)
    _check_process_bounds(get_estimates(estimate), p, q)
    nothing
end

function _check_process_bounds(estimate::AbstractArray{<:Number, M}, p, q) where {M}
    checkbounds(Bool, axes(estimate, 1), p) ||
        throw(BoundsError(
            "The first process index, $p, is out of bounds for estimate with " *
            "size $(size(estimate))"
        ))
    checkbounds(Bool, axes(estimate, 2), q) ||
        throw(BoundsError(
            "The second process index, $q, is out of bounds for estimate with " *
            "size $(size(estimate))"
        ))
end

function _check_process_bounds(
        ::AbstractArray{SMatrix{P, Q, T, L}, D}, p, q) where {P, Q, T, L, D}
    checkbounds(Bool, 1:P, p) ||
        throw(BoundsError(
            "The first process index, $p, is out of bounds for estimate with " *
            "size $((P, Q))"
        ))
    checkbounds(Bool, 1:Q, q) ||
        throw(BoundsError(
            "The second process index, $q, is out of bounds for estimate with " *
            "size $((P, Q))"
        ))
end

function _check_process_bounds(
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}}, p, q) where {N}
    # If one is an Int this still iterates, if both are Int then method below is used
    [_check_process_bounds(estimate, pᵢ, qᵢ) for (pᵢ, qᵢ) in Iterators.product(p, q)]
end

function _check_process_bounds(
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}}, p::Int, q::Int) where {N}
    (min(p, q), max(p, q)) in keys(estimate) ||
        throw(BoundsError(
            "Process indices ($p, $q) are out of bounds for estimate with " *
            "keys $(keys(estimate))"
        ))
    nothing
end

# Index bounds checking
function check_index_bounds(
        estimate::AbstractEstimate{E, D, N}, i::Vararg{Any, N}) where {E, D, N}
    _check_index_bounds(get_estimates(estimate), i...)
end

function _check_index_bounds(
        estimate::AbstractArray{<:Number, M}, i::Vararg{Any, N}) where {M, N}
    checkbounds(
        Bool,
        selectdim(
            selectdim(estimate, 1, first(axes(estimate, 1))),
            1,
            first(axes(estimate, 2))
        ),
        i...
    ) || throw(BoundsError(
        "Index $i is out of bounds for estimate with size $(size(estimate)[3:end])"
    ))
end

function _check_index_bounds(estimate::AbstractArray{SMatrix{P, Q, T, L}, D},
        i::Vararg{Any, D}) where {P, Q, T, L, D}
    checkbounds(Bool, estimate, i...) ||
        throw(BoundsError(
            "Index $i is out of bounds for estimate with size $(size(estimate))"
        ))
end

function _check_index_bounds(
        estimate::Dict{<:Tuple, <:AbstractArray{<:Number, N}},
        i::Vararg{Any, N}) where {N}
    if !all(x -> checkbounds(Bool, x, i...), values(estimate))
        all(size(x) == size(first(values(estimate))), values(estimate)) ||
            throw(BoundsError(
                "In an estimate, all estimates in the dictionary should have " *
                "the same size, these do not."
            ))
        throw(BoundsError(
            "Index $i is out of bounds for estimate with " *
            "size $(size(first(values(estimate))))"
        ))
    end
    nothing
end

# Getindex implementation
"""
    Base.getindex(estimate::AbstractEstimate, p, q)
    Base.getindex(estimate::AbstractEstimate, p, q, i...)

Get a subset of the estimate corresponding to processes `p` and `q` and
potentially spatial indices `i`.

If `p` and `q` are integers, the result will be a single process estimate.
Otherwise the result will be the same type of estimate but with the specified
processes. If one of the indices is an integer and the other is not, the result
will still have the original number of dimensions internally.
"""
function Base.getindex(
        estimate::AbstractEstimate{E, D, N}, p, q, i::Vararg{Any, N}) where {E, D, N}
    checkbounds(estimate, p, q, i...)
    return _construct_estimate_subset(
        typeof(estimate),
        E,
        get_evaluation_points_index(estimate, i...),
        get_estimates_index(estimate, p, q, i...),
        get_process_information_index(estimate, p, q),
        get_estimation_information(estimate),
        get_extra_information(estimate)...
    )
end

function Base.getindex(estimate::AbstractEstimate{E}, p, q) where {E}
    checkbounds(estimate, p, q)
    return _construct_estimate_subset(
        typeof(estimate),
        E,
        get_evaluation_points(estimate),
        get_estimates_index(estimate, p, q),
        get_process_information_index(estimate, p, q),
        get_estimation_information(estimate),
        get_extra_information(estimate)...
    )
end

# Fallback for types that need trait parameter explicitly
function _construct_estimate_subset(
        ::Type{T}, trait::Type{<:EstimateTrait}, args...) where {T <: AbstractEstimate}
    # Extract the base constructor and call it with the trait
    base_constructor = constructorof(T)
    return base_constructor{trait}(args...)
end

# Helper functions for indexing
function get_evaluation_points_index(
        estimate::AbstractEstimate{E, D, N}, i::Vararg{Any, N}) where {E, D, N}
    _get_argument_index(get_evaluation_points(estimate), i...)
end

function _get_argument_index(argument::NTuple{N}, i::Vararg{Any, N}) where {N}
    getindex.(argument, tuple(i...))
end

function _get_argument_index(argument, i)
    argument[i]
end

function get_estimates_index(
        estimate::AbstractEstimate{E, D, N}, p, q, i::Vararg{Any, N}) where {E, D, N}
    return _get_estimate_index(process_trait(estimate), get_estimates(estimate), p, q, i...)
end

function get_estimates_index(estimate::AbstractEstimate{E, D, N}, p, q) where {E, D, N}
    return _get_estimate_index(process_trait(estimate), get_estimates(estimate), p, q)
end

# Estimate indexing implementations for different process traits
function _get_estimate_index(::SingleProcessTrait, estimate::AbstractArray{<:Number, N},
        p::Int, q::Int, i::Vararg{Any, N}) where {N}
    return estimate[i...]
end

function _get_estimate_index(::SingleProcessTrait, estimate::AbstractArray{<:Number},
        p::Int, q::Int)
    return estimate
end

function _get_estimate_index(::MultipleTupleTrait, estimate::AbstractArray{<:SMatrix, N},
        p, q, i::Vararg{Any, N}) where {N}
    return _get_estimate_index(MultipleTupleTrait(), estimate, p, q)[i...]
end

function _get_estimate_index(::MultipleTupleTrait, estimate::AbstractArray{<:SMatrix, N},
        p, q) where {N}
    return getindex.(estimate, Ref(_index_to_svec(p)), Ref(_index_to_svec(q)))
end

function _get_estimate_index(::MultipleTupleTrait, estimate::AbstractArray{<:SMatrix, N},
        p::Int, q::Int) where {N}
    return getindex.(estimate, p, q)
end

_index_to_svec(i::Int) = SVector(i)
_index_to_svec(i::SVector) = i
_index_to_svec(i::SOneTo) = i
_index_to_svec(i::AbstractVector) = SVector(i...)

function _get_estimate_index(::MultipleVectorTrait, estimate::AbstractArray{<:Number, M},
        p::Int, q::Int, i::Vararg{Any, N}) where {M, N}
    estimate[p, q, i...]
end

function _get_estimate_index(::MultipleVectorTrait, estimate::AbstractArray{<:Number, M},
        p, q, i::Vararg{Any, N}) where {M, N}
    estimate[_index_to_vec(p), _index_to_vec(q), i...]
end

_index_to_vec(i::Int) = i:i
_index_to_vec(i::AbstractVector) = i

function _get_estimate_index(::MultipleVectorTrait, estimate::AbstractArray{<:Number, N},
        p::Int, q::Int) where {N}
    # Once first is selected, second becomes first
    estimate[p, q, ntuple(Returns(:), Val{N - 2}())]
end

function _get_estimate_index(::MultipleVectorTrait, estimate::AbstractArray{<:Number, N},
        p, q) where {N}
    # Once first is selected, second becomes first
    estimate[_index_to_vec(p), _index_to_vec(q), ntuple(Returns(:), Val{N - 2}())...]
end
