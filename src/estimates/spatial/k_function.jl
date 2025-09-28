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

k_function(data, region; kwargs...) = k_function(spatial_data(data, region); kwargs...)
function k_function(data::SpatialData; kwargs...)
    return k_function(c_function(data; kwargs...))
end
function k_function(c::CFunction{E, D}) where {E, D}
    mean_prod = getprocessinformation(c).mean_product
    radii = getargument(c)
    value = C2K(getargument(c), getestimate(c), process_trait(c), mean_prod, Val{D}())
    processinfo = getprocessinformation(c)
    estimationinfo = getestimationinformation(c)
    return KFunction{E}(radii, value, processinfo, estimationinfo)
end
function k_function(spectrum::NormalOrRotationalSpectra; kwargs...)
    k_function(c_function(spectrum; kwargs...))
end

function partial_k_function(data, region; kwargs...)
    partial_k_function(spatial_data(data, region); kwargs...)
end
function partial_k_function(data::SpatialData; kwargs...)
    k_function(partial_c_function(data; kwargs...))
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
function C2K(radii, c::AbstractArray,
        ::MultipleVectorTrait, mean_prod, ::Val{D}) where {D}
    out = zeros(eltype(c), size(c))
    for idx in CartesianIndices(size(c)[1:(ndims(c) - 1)])
        mean_prod_slice = mean_prod[idx]
        for (i, radius) in enumerate(radii)
            out[idx, i] = _C2K(radius, c[idx, i], mean_prod_slice, Val{D}())
        end
    end
    return out
end

function C2K(radii, c::AbstractArray,
        ::Union{MultipleTupleTrait, SingleProcessTrait}, mean_prod, ::Val{D}) where {D}
    out = zeros(eltype(c), size(c))
    for (i, radius) in enumerate(radii)
        out[i] = _C2K(radius, c[i], mean_prod, Val{D}())
    end
    return out
end

function _C2K(radius, c, mean_prod, ::Val{D}) where {D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    return c ./ mean_prod .+ (V * (radius^D))
end

function _C2K(radius, c, mean_prod, ::Val{2})
    return c ./ mean_prod .+ (pi * radius^2)
end
