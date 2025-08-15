function check_spectral_indices(
    indices,
    input::Union{
        SpectralEstimate{D,F,P,N},
        PartialSpectralEstimate{D,F,P,N},
        NTuple{D,Union{GeoTable,PointSet}},
    },
) where {D,F,P,N}
    @assert indices isa
            AbstractVector{<:Tuple{Int,Int,AbstractVector{Int},AbstractVector{Int}}} ||
            indices isa AbstractVector{Tuple{Int,Int}}
    for index in indices
        for i = 1:D
            @assert index[i] ⊆ 1:D "problem with the indices, should all be subsets of $(1:D), but got $(index)"
        end
    end
end
