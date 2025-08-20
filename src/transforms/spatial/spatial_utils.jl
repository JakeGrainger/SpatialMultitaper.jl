# TODO: this will be redundent soon
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
    indices::Nothing,
    input::Union{
        SpectralEstimate{D,F,P,N},
        PartialSpectra{D,F,P,N},
        NTuple{P,Union{GeoTable,PointSet}},
    },
) where {D,F,P,N}
    nothing
end

function check_spectral_indices(indices::Nothing, input::Union{GeoTable,PointSet})
    nothing
end


function default_indices(
    input::Union{
        SpectralEstimate{D,F,P,AbstractArray{<:Number,M}},
        PartialSpectra{D,F,P,AbstractArray{<:Number,M}},
    },
) where {D,F,P,M}
    return [(i, j) for i = 1:P for j = i:P]
end

function default_indices(
    input::Union{SpectralEstimate{D,F,P,N},PartialSpectra{D,F,P,N}},
) where {D,F,P,N}
    return nothing
end

function default_indices(input::AbstractVector{Union{GeoTable,PointSet}})
    return [(i, j) for i = 1:length(input) for j = i:length(input)]
end

function default_indices(input::NTuple{P,Union{GeoTable,PointSet}}) where {P}
    return nothing
end

function default_indices(
    input::Union{
        SpectralEstimate{D,F,1,AbstractArray{<:Number,M}},
        PartialSpectra{D,F,1,AbstractArray{<:Number,M}},
        GeoTable,
        PointSet,
    },
) where {D,F,M}
    return nothing
end