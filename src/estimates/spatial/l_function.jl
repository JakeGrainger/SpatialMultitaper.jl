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

l_function(data, region; kwargs...) = l_function(spatial_data(data, region); kwargs...)
function l_function(data::SpatialData; kwargs...)
    return l_function(k_function(data; kwargs...))
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
    partial_l_function(spatial_data(data, region); kwargs...)
end
function partial_l_function(data::SpatialData; kwargs...)
    l_function(partial_k_function(data; kwargs...))
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
function K2L(k::AbstractArray, ::Val{D}) where {D}
    out = zeros(eltype(k), size(k))
    for idx in eachindex(out)
        out[idx] = _K2L(k[idx], Val{D}())
    end
    return out
end

function _K2L(k, ::Val{D}) where {D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    sign.(k) .* (abs.(k) ./ V) .^ (1 / D) # broadcast because may contain SVectors
end

function _K2L(k, ::Val{2})
    sign.(k) .* sqrt.(abs.(k) ./ pi) # broadcast because may contain SVectors
end
