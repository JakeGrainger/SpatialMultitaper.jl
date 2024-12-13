abstract type FrequencyDomainEstimate{P} end
# should be a struct implementing getfreq and getestimate

function Base.getindex(est::FrequencyDomainEstimate{1}, i, j)
    i == j == 1 || throw(BoundsError("Estimate is scalar, only index (1, 1) is valid"))
    return est
end

function Base.getindex(est::FrequencyDomainEstimate{P}, i::Int, j::Int) where {P}
    1 <= i <= P || throw(BoundsError("Index $i is out of bounds for estimate with $P processes"))
    1 <= j <= P || throw(BoundsError("Index $j is out of bounds for estimate with $P processes"))
    return constructorof(typeof(est))(getfreq(est), getestindex(getestimate(est), i, j))
end

function getestindex(est::AbstractArray{S,D}, i, j) where {S<:StaticMatrix, D}
    getindex.(est, i, j)
end

function getestindex(est::AbstractArray{T,D}, i, j) where {T<:Number, D}
    est_flat = reshape(est, size(est)[1:2]..., prod(size(est)[3:end]))
    est_subset = est_flat[i, j, :]
    return reshape(est_subset, size(est)[3:end])
end

function checkfreqdomaininputs(
        freq::NTuple{D, F},
		est::AbstractArray{SMatrix{P, P, T, L}, N},
	) where {D, F, P, T, L, N}
        @assert N == D "frequencies and estimate should be same dimension when estimate is an array of matrices."
		@assert length.(freq) == size(est)
		return P
end

function checkfreqdomaininputs(
    freq::NTuple{D, F},
    est::AbstractArray{T, D},
) where {D, F, T}
    @assert length.(freq) == size(est)
    return 1
end

function checkfreqdomaininputs(
    freq::NTuple{D, F},
    est::AbstractArray{T, N},
) where {D, F, T, N}
    @assert N == D+2 "frequencies should either have same dimension as estimate (if one process being analysed) or the dim(estimate)==dim(freq)+2."
    @assert length.(freq) == size(est)[3:end]
    return size(est, 1)
end
