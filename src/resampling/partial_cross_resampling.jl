struct PartialCrossResampler{D, S}
    data::D
    shift_method::S
end

function Base.rand(rng::AbstractRNG, resampler::PartialCrossResampler)
    return shifted_data(rng, resampler.data, resampler.shift_method)
end

function shifted_data(rng, data::NTuple{P}, shift_method) where {P}
    ntuple(p -> marginal_shift(data[p], rand(rng, shift_method)), Val{P - 1}()) # dont need last one as use symmetry
end
