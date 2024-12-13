abstract type FrequencyDomainEstimate{D,P} end
# should be a struct implementing getfreq and getestimate

function checkprocessbounds(::FrequencyDomainEstimate{D,P}, p) where {D,P}
    checkbounds(Bool, 1:P, p) ||
        throw(BoundsError("Index $p is out of bounds for estimate with $P processes"))
    nothing
end
function Base.checkbounds(
    est::FrequencyDomainEstimate{D,P},
    p,
    q,
    i::Vararg{Any,D},
) where {D,P}
    checkprocessbounds(est, p)
    checkprocessbounds(est, q)
    # TODO should check i and have option to disable
    nothing
end
function Base.checkbounds(est::FrequencyDomainEstimate, p, q)
    checkprocessbounds(est, p)
    checkprocessbounds(est, q)
    nothing
end

function Base.getindex(
    est::FrequencyDomainEstimate{D,P},
    p,
    q,
    i::Vararg{Any,D},
) where {D,P}
    checkbounds(est, p, q, i...)
    return constructorof(typeof(est))(getfreqindex(est, i...), getestimateindex(est, p, q, i...))
end

function Base.getindex(est::FrequencyDomainEstimate{D,P}, p::Int, q::Int) where {D,P}
    checkbounds(est, p, q)
    return constructorof(typeof(est))(getfreq(est), getestimateindex(est, p, q))
end

function getfreqindex(est::FrequencyDomainEstimate{D,P}, i::Vararg{Any,D}) where {D,P}
    getindex.(getfreq(est), i)
end

function getestimateindex(est::FrequencyDomainEstimate, p, q, i...)
    _getestimateindex(getestimate(est), p, q, i...)
end

function _getestimateindex(
    est::AbstractArray{<:SMatrix,D},
    p,
    q,
    i::Vararg{Any,D},
) where {D}
    _getestimateindex(est, p, q)[i...]
end
function _getestimateindex(est::AbstractArray{<:SMatrix,D}, p, q) where {D}
    getindex.(est, p, q)
end

function _getestimateindex(
    est::AbstractArray{<:Number,N},
    p,
    q,
    i::Vararg{Any,D},
) where {N,D}
    est[p, q, i...]
end

function _getestimateindex(est::AbstractArray{<:Number,N}, p, q) where {N}
    collect(selectdim(est(selectdim(est, 1, p)), 2, q))
end

function _getestimateindex(est::AbstractArray{<:Number,D}, p, q, i::Vararg{Any,D}) where {D} # one process case, p and q have been checked by this point but should be 1
    est[i...]
end


function checkfreqdomaininputs(
    freq::NTuple{D,F},
    est::AbstractArray{SMatrix{P,P,T,L},N},
) where {D,F,P,T,L,N}
    @assert N == D "frequencies and estimate should be same dimension when estimate is an array of matrices."
    @assert length.(freq) == size(est)
    return P
end

function checkfreqdomaininputs(freq::NTuple{D,F}, est::AbstractArray{T,D}) where {D,F,T}
    @assert length.(freq) == size(est)
    return 1
end

function checkfreqdomaininputs(freq::NTuple{D,F}, est::AbstractArray{T,N}) where {D,F,T,N}
    @assert N == D + 2 "frequencies should either have same dimension as estimate (if one process being analysed) or the dim(estimate)==dim(freq)+2."
    @assert length.(freq) == size(est)[3:end]
    return size(est, 1)
end
