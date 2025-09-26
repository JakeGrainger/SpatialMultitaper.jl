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
    value = L2centeredL(getargument(l), getestimate(l))
    processinfo = getprocessinformation(k)
    estimationinfo = getestimationinformation(k)
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

# TODO: only works for array of SMatrix
function L2centeredL(radii, k_function)
    [k .- r for (k, r) in zip(k_function, radii)]
end
