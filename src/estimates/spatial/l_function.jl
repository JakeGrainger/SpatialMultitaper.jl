struct LFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function LFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        P, Q = checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        new{E, D, P, Q, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end
getbaseestimatename(::Type{<:LFunction}) = "L function"
getargument(f::LFunction) = f.radii
getestimate(f::LFunction) = f.value

function l_function(data, region; kwargs...)
    return l_function(k_function(data, region; kwargs...))
end

l_function(c::CFunction) = l_function(k_function(c))
function l_function(spec::NormalOrRotationalSpectra; kwargs...)
    l_function(k_function(spec; kwargs...))
end

function l_function(k::KFunction{E, D}) where {E, D}
    radii = getargument(k)
    value = K2L(getestimate(k), Val{D}())
    processinfo = getprocessinformation(k)
    estimationinfo = getestimationinformation(k)
    return LFunction{E}(radii, value, processinfo, estimationinfo)
end

function partial_l_function(data, region; kwargs...)
    l_function(partial_k_function(data, region; kwargs...))
end
function partial_l_function(spectrum::NormalOrRotationalSpectra{PartialTrait}; kwargs...)
    k_function(spectrum; kwargs...)
end
function partial_l_function(spectrum::NormalOrRotationalSpectra{MarginalTrait}; kwargs...)
    k_function(partial_spectra(spectrum); kwargs...)
end

partial_l_function(est::CFunction{PartialTrait}) = l_function(est)
partial_l_function(est::KFunction{PartialTrait}) = l_function(est)
function partial_l_function(est::Union{CFunction{MarginalTrait}, KFunction{MarginalTrait}})
    throw(partial_from_marginal_error(LFunction, typeof(est)))
end

## internals

function K2L(k_function, ::Val{D}) where {D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    return [sign.(k) .* (abs.(k) ./ V) .^ (1 / D) for k in k_function] # broadcast because may contain SVectors
end

function K2L(k_function, ::Val{2})
    return [sign.(k) .* sqrt.(abs.(k) ./ pi) for k in k_function] # broadcast because may contain SVectors
end

function L_function(k::KFunction{R, K, I, T, D, P, Q}) where {R, K, I, T, D, P, Q}
    return LFunction(k.radii, K2L(k.K_function, Val{D}()),
        getprocessinformation(k), getestimationinformation(k))
end

function L_function(
        data,
        region;
        radii,
        nfreq,
        fmax,
        tapers,
        mean_method::MeanEstimationMethod = DefaultMean(),
        freq_radii = default_rotational_radii(nfreq, fmax),
        rotational_method = default_rotational_kernel(freq_radii)
)
    k = K_function(
        data,
        region;
        radii = radii,
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
        freq_radii = freq_radii,
        rotational_method = rotational_method
    )
    return L_function(k)
end
