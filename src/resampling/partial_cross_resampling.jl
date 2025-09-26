struct PartialCrossResampler{D, S}
    data::D
    shift_method::S
end

function Base.rand(rng::AbstractRNG, resampler::PartialCrossResampler)
    return shifted_data(rng, resampler.data, resampler.shift_method)
end

function shifted_data(rng, data::MultipleSpatialDataTuple{P}, shift_method) where {P}
    spatial_data(
        ntuple(p -> observations(marginal_shift(data[p], rand(rng, shift_method))),
            Val{P - 1}()),
        getregion(data))
end
