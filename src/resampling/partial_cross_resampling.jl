struct PartialCrossResampler{D, R, S, K}
    data::D
    region::R
    shift_method::S
    mt_kwargs::K
end

function Base.rand(rng::AbstractRNG, resampler::PartialCrossResampler)
    region = resampler.region
    tapers = resampler.mt_kwargs.tapers
    nfreq = resampler.mt_kwargs.nfreq
    fmax = resampler.mt_kwargs.fmax
    mean_method = resampler.mt_kwargs.mean_method

    new_data = shifted_data(rng, resampler.data, resampler.shift_method)

    return tapered_dft(new_data, tapers, nfreq, fmax, region, mean_method)
end

function shifted_data(rng, data::NTuple{P}, shift_method) where {P}
    ntuple(p -> marginal_shift(data[p], rand(rng, shift_method)), Val{P}())
end
