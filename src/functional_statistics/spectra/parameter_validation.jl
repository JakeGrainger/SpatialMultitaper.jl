#### single parameter validation functions
function validate_dk(dk)
    @argcheck all(@. dk isa Number)
    @argcheck all(@. dk > 0)
end
function validate_nk(nk)
    @argcheck all(@. nk > 0)
end
function validate_kmax(kmax)
    @argcheck all(@. kmax isa Number)
    @argcheck all(@. kmax > 0)
end
function validate_nw(nw)
    @argcheck all(@. nw isa Number)
    @argcheck all(@. nw > 0)
end

function validate_tapers(tapers)
    nothing
end

function validate_mean_estimation_method(mean_method)
    @argcheck mean_method isa MeanEstimationMethod
end

#### resolution of dk/nk/kmax

function resolve_dk_nk_kmax(data::SpatialData; kwargs...)
    nk = get(kwargs, :nk, nothing)
    kmax = get(kwargs, :kmax, nothing)
    dk = get(kwargs, :dk, nothing)
    _nk, _kmax = _resolve_dk_nk_kmax(nk, kmax, dk, data)

    other_kwargs = filter(kv -> !(kv[1] in (:nk, :kmax, :dk)), kwargs)
    return (; nk = _nk, kmax = _kmax, other_kwargs...)
end

function _resolve_dk_nk_kmax(nk, kmax, ::Nothing, ::SpatialData{T, D}) where {T, D}
    _nk = _resolve_nk(nk, Val{D}())
    _kmax = _resolve_kmax(kmax, Val{D}())
    return _nk, _kmax
end

function _resolve_dk_nk_kmax(nk, kmax, dk, data::SpatialData{T, D}) where {T, D}
    @warn "Specifying all three of nk, kmax, and dk is not necessary. Only two are needed. Ignoring dk."
    return _resolve_dk_nk_kmax(nk, kmax, nothing, data)
end

function _resolve_dk_nk_kmax(nk, ::Nothing, dk, ::SpatialData{T, D}) where {T, D}
    _nk = _resolve_nk(nk, Val{D}())
    _dk = _resolve_dk(dk, Val{D}())
    _kmax = _nk .* _dk ./ 2
    return _nk, _kmax
end

function _resolve_dk_nk_kmax(::Nothing, kmax, dk, ::SpatialData{T, D}) where {T, D}
    _kmax = _resolve_kmax(kmax, Val{D}())
    _dk = _resolve_dk(dk, Val{D}())
    _nk = ceil.(Int, 2 .* _kmax ./ _dk)
    _nk = _resolve_nk(_nk, Val{D}())
    return _nk, _kmax
end

function _resolve_dk_nk_kmax(nk::Nothing, kmax, ::Nothing, data::SpatialData)
    dk = default_dk(data, nk, kmax)
    return _resolve_dk_nk_kmax(nk, kmax, dk, data)
end

function _resolve_dk_nk_kmax(nk, kmax::Nothing, ::Nothing, data::SpatialData)
    dk = default_dk(data, nk, kmax)
    return _resolve_dk_nk_kmax(nk, kmax, dk, data)
end

function _resolve_dk_nk_kmax(
        nk::Nothing, kmax::Nothing, dk::Nothing, data::SpatialData)
    error("Must specify at least two of nk, kmax, or dk.")
end

function _resolve_dk_nk_kmax(nk::Nothing, kmax::Nothing, dk, data::SpatialData)
    error("Must specify at least two of nk, kmax, or dk, not just dk.")
end

function _resolve_nk(nk, ::Val{D}) where {D}
    _nk = floor.(Integer, nk)
    if !all(@. nk == _nk)
        @warn "Non-integer nk rounded down to nearest integer."
    end
    if D > 1 && any(@. _nk >= 1000)
        @warn "Number of wavenumbers `nk` >= 1000, probably means you misspecified kmax, may want to ctrl-c and check."
    end
    @argcheck all(@. _nk > 0)
    return _resolve_dims(_nk, Val{D}())
end

function _resolve_dk(dk, ::Val{D}) where {D}
    @argcheck all(@. dk isa Number)
    @argcheck all(@. dk > 0)
    return _resolve_dims(dk, Val{D}())
end

function _resolve_kmax(kmax, ::Val{D}) where {D}
    @argcheck all(@. kmax isa Number)
    @argcheck all(@. kmax > 0)
    return _resolve_dims(kmax, Val{D}())
end

function _resolve_dims(x::AbstractVector, ::Val{D}) where {D}
    return _resolve_dims(tuple(x...), Val{D}())
end
function _resolve_dims(x::NTuple{D}, ::Val{D}) where {D}
    return x
end
function _resolve_dims(x::NTuple{L}, ::Val{D}) where {D, L}
    throw(ArgumentError("Expected a tuple of length $D, got length $L."))
end
function _resolve_dims(x::Number, ::Val{D}) where {D}
    return ntuple(Returns(x), Val{D}())
end

default_dk(data::SpatialData, nk, kmax) = default_dk(getregion(data), nk, kmax)
default_dk(region::Geometry, ::Nothing, kmax) = _default_dk(region)
default_dk(region::Geometry, nk, ::Nothing) = _default_dk(region)
default_dk(region::Geometry, ::Nothing, ::Nothing) = _default_dk(region)
default_dk(region::Geometry, nk, kmax) = nothing
function _default_dk(region::Geometry)
    bbox = boundingbox(region)
    return 1 ./ Meshes.ustrip.(sides(bbox))
end

#### taper resolution
function resolve_tapers(data; kwargs...)
    region = getregion(data)
    tapers = get(kwargs, :tapers, nothing)
    nw = get(kwargs, :nw, nothing)
    tapers = _resolve_tapers(region, tapers, nw)
    other_kwargs = filter(kv -> !(kv[1] in (:tapers, :nw)), kwargs)
    return (; tapers = tapers, other_kwargs...)
end
function _resolve_tapers(region, tapers::TaperFamily, nw::Nothing)
    return tapers
end
function _resolve_tapers(region, tapers::TaperFamily, nw)
    @warn "Both `tapers` and `nw` specified, ignoring `nw`."
    return tapers
end
function _resolve_tapers(region, ::Nothing, ::Nothing)
    nw = 4
    return _resolve_tapers(region, nothing, nw)
end
function _resolve_tapers(region, ::Nothing, nw)
    return _default_tapers(region, nw)
end

function _default_tapers(region::Box, nw)
    M = floor(Int, 2 .* nw .- 1)
    M = _resolve_dims(M, Val{embeddim(region)}())
    return sin_taper_family(M, region)
end
function _default_tapers(region, nw)
    bandwidth = nw / Meshes.ustrip(minimum(sides(boundingbox(region))))
    return make_tapers(region, bandwidth = bandwidth)
end

function resolve_mean_method(data; kwargs...)
    mean_method = get(kwargs, :mean_method, nothing)
    mean_method = _resolve_mean_method(mean_method)
    other_kwargs = filter(kv -> kv[1] != :mean_method, kwargs)
    return (; mean_method = mean_method, other_kwargs...)
end
function _resolve_mean_method(::Nothing)
    return DefaultMean()
end
function _resolve_mean_method(mean_method)
    return mean_method
end
