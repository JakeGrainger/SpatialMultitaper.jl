struct CenteredLFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function CenteredLFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        P, Q = checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        new{E, D, P, Q, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end
getbaseestimatename(::Type{<:CenteredLFunction}) = "centered L function (L(r)-r)"
getargument(f::CenteredLFunction) = f.radii
getestimate(f::CenteredLFunction) = f.value

function centered_l_function(data, region; kwargs...)
    centered_l_function(spatial_data(data, region); kwargs...)
end
function centered_l_function(data::SpatialData; kwargs...)
    return centered_l_function(l_function(data; kwargs...))
end
function centered_l_function(spectrum::NormalOrRotationalSpectra; kwargs...)
    centered_l_function(l_function(spectrum; kwargs...))
end
centered_l_function(c::CFunction) = centered_l_function(l_function(c))
centered_l_function(k::KFunction) = centered_l_function(l_function(k))
function centered_l_function(l::LFunction{E, D}) where {E, D}
    radii = getargument(l)
    value = L2centeredL(getargument(l), getestimate(l), process_trait(l))
    processinfo = getprocessinformation(l)
    estimationinfo = getestimationinformation(l)
    return CenteredLFunction{E}(radii, value, processinfo, estimationinfo)
end

function partial_centered_l_function(data, region; kwargs...)
    partial_centered_l_function(spatial_data(data, region); kwargs...)
end
function partial_centered_l_function(data::SpatialData; kwargs...)
    centered_l_function(partial_l_function(data; kwargs...))
end
function partial_centered_l_function(
        spectrum::NormalOrRotationalSpectra{PartialTrait}; kwargs...)
    k_function(spectrum; kwargs...)
end
function partial_centered_l_function(
        spectrum::NormalOrRotationalSpectra{MarginalTrait}; kwargs...)
    k_function(partial_spectra(spectrum); kwargs...)
end
partial_centered_l_function(est::CFunction{PartialTrait}) = centered_l_function(est)
partial_centered_l_function(est::KFunction{PartialTrait}) = centered_l_function(est)
partial_centered_l_function(est::LFunction{PartialTrait}) = centered_l_function(est)
function partial_centered_l_function(est::Union{
        LFunction{MarginalTrait}, CFunction{MarginalTrait}, KFunction{MarginalTrait}})
    throw(partial_from_marginal_error(LFunction, typeof(est)))
end

## internals

function L2centeredL(radii, est::AbstractArray, ::MultipleVectorTrait)
    out = zeros(eltype(est), size(est))
    for idx in CartesianIndices(size(est)[1:(ndims(est) - 1)])
        for (i, radius) in enumerate(radii)
            out[idx, i] = _L2centeredL(radius, est[idx, i])
        end
    end
    return out
end

function L2centeredL(
        radii, est::AbstractArray, ::Union{MultipleTupleTrait, SingleProcessTrait})
    out = zeros(eltype(est), size(est))
    for (i, radius) in enumerate(radii)
        out[i] = _L2centeredL(radius, est[i])
    end
    return out
end

function _L2centeredL(radius, est)
    est .- radius
end
