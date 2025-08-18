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

function check_spectral_indices(
    indices,
    input::Union{
        SpectralEstimate{D,F,1,N},
        PartialSpectra{D,F,1,N},
        NTuple{1,Union{GeoTable,PointSet}},
        GeoTable,
        PointSet,
    },
) where {D,F,N}
    @assert indices isa Nothing "in the univariate case, indices should be `nothing`"
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

function default_indices(
    input::Union{
        SpectralEstimate{D,F,1,N},
        PartialSpectra{D,F,1,N},
        NTuple{1,Union{GeoTable,PointSet}},
        GeoTable,
        PointSet,
    },
) where {D,F,N}
    return nothing
end