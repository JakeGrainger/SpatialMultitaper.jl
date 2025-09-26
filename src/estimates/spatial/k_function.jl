struct KFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function KFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        P, Q = checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        new{E, D, P, Q, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end
getbaseestimatename(::Type{<:KFunction}) = "K function"
getargument(f::KFunction) = f.radii
getestimate(f::KFunction) = f.value

function k_function(data, region; kwargs...)
    return k_function(c_function(data, region; kwargs...))
end

function k_function(c::CFunction{E, D}) where {E, D}
    mean_prod = getprocessinformation(c).mean_product
    radii = getargument(c)
    value = C2K(getargument(c), getestimate(c), mean_prod, Val{D}())
    processinfo = getprocessinformation(c)
    estimationinfo = getestimationinformation(c)
    return KFunction{E}(radii, value, processinfo, estimationinfo)
end
function k_function(spectrum::NormalOrRotationalSpectra; kwargs...)
    k_function(c_function(spectrum; kwargs...))
end

function partial_k_function(data, region; kwargs...)
    k_function(partial_c_function(data, region; kwargs...))
end
function partial_k_function(spectrum::NormalOrRotationalSpectra{PartialTrait}; kwargs...)
    k_function(spectrum; kwargs...)
end
function partial_k_function(spectrum::NormalOrRotationalSpectra{MarginalTrait}; kwargs...)
    k_function(partial_spectra(spectrum); kwargs...)
end
partial_k_function(c::CFunction{PartialTrait}) = k_function(c)
function partial_k_function(::CFunction{MarginalTrait})
    throw(partial_from_marginal_error(KFunction, CFunction))
end

# internals

function C2K(radii, c_function, mean_prod, ::Val{D}) where {D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    return [c ./ mean_prod .+ (V * (r^D)) for (r, c) in zip(radii, c_function)]
end

function C2K(radii, c_function, mean_prod, ::Val{2})
    return [c ./ mean_prod .+ (pi * r^2) for (r, c) in zip(radii, c_function)]
end
