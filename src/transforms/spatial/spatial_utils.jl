function check_spectral_indices(
    indices,
    input::Union{
        SpectralEstimate{D,F,P,N},
        PartialSpectra{D,F,P,N},
        NTuple{P,Union{GeoTable,PointSet}},
    },
) where {D,F,P,N}
    @assert indices isa
            AbstractVector{<:Tuple{Int,Int,AbstractVector{Int},AbstractVector{Int}}} ||
            indices isa AbstractVector{Tuple{Int,Int}}
    for index in indices
        for i in index
            @assert i ⊆ 1:P "problem with the indices, should all be subsets of $(1:P), but got $(index)"
        end
    end
end

function default_indices(
    input::Union{
        SpectralEstimate{D,F,P,N},
        PartialSpectra{D,F,P,N},
        NTuple{P,Union{GeoTable,PointSet}},
    },
) where {D,F,P,N}
    return [(i, j) for i = 1:P for j = i:P]
end