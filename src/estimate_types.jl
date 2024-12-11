abstract type FrequencyDomainEstimate{P} end
# should be a struct implementing getfreq and getestimate

function Base.getindex(est::FrequencyDomainEstimate{1}, i, j)
    i == j == 1 || throw(BoundsError("Estimate is scalar, only index (1, 1) is valid"))
    return est
end

function Base.getindex(est::FrequencyDomainEstimate{P}, i::Int, j::Int) where {P}
    1 <= i <= P || throw(BoundsError("Index $i is out of bounds for estimate with $P processes"))
    1 <= j <= P || throw(BoundsError("Index $j is out of bounds for estimate with $P processes"))
    return constructorof(est)(getfreq(est), getindex.(getestimate(est), i, j))
end

function checkfreqdomaininputs(
        freq,
		est::AbstractArray{SMatrix{P, P, T, L}, D},
	) where {P, T, L, D}
		@assert length.(freq) == size(est)
		return P
end

function checkfreqdomaininputs(
    freq,
    est::AbstractArray{T, D},
) where {T, D}
    @assert length.(freq) == size(est)
    return 1
end
