# defaults
default_dk(data::SpatialData, nk, kmax) = default_dk(getregion(data), nk, kmax)
default_dk(region::Geometry, ::Nothing, kmax) = _default_dk(region)
default_dk(region::Geometry, nk, ::Nothing) = _default_dk(region)
default_dk(region::Geometry, ::Nothing, ::Nothing) = _default_dk(region)
default_dk(region::Geometry, nk, kmax) = nothing
function _default_dk(region::Geometry)
    bbox = boundingbox(region)
    return 1 ./ Meshes.ustrip.(sides(bbox))
end

# validation
function _validate_dims(x::AbstractVector, ::Val{D}) where {D}
    return _validate_dims(tuple(x...), Val{D}())
end
function _validate_dims(x::NTuple{D}, ::Val{D}) where {D}
    return x
end
function _validate_dims(x::NTuple{L}, ::Val{D}) where {D, L}
    throw(ArgumentError("Expected a tuple of length $D, got length $L."))
end
function _validate_dims(x::Number, ::Val{D}) where {D}
    return ntuple(Returns(x), Val{D}())
end

function _validate_nk(nk, ::Val{D}) where {D}
    _nk = floor.(Integer, nk)
    if !all(@. nk == _nk)
        @warn "Non-integer nk rounded down to nearest integer."
    end
    if D > 1 && any(@. _nk >= 1_000)
        @warn "Number of wavenumbers `nk` >= 1000, probably means you misspecified kmax, may want to ctrl-c and check."
    end
    @argcheck all(@. _nk > 0)
    return _validate_dims(_nk, Val{D}())
end

function _validate_dk(dk, ::Val{D}) where {D}
    @argcheck all(@. dk isa Number)
    @argcheck all(@. dk > 0)
    return _validate_dims(dk, Val{D}())
end

function _validate_kmax(kmax, ::Val{D}) where {D}
    @argcheck all(@. kmax isa Number)
    @argcheck all(@. kmax > 0)
    return _validate_dims(kmax, Val{D}())
end

function _validate_wavenumber_params(nk, kmax, ::Nothing, ::SpatialData{T, D}) where {T, D}
    _nk = _validate_nk(nk, Val{D}())
    _kmax = _validate_kmax(kmax, Val{D}())
    return _nk, _kmax
end
function _validate_wavenumber_params(nk, ::Nothing, dk, ::SpatialData{T, D}) where {T, D}
    _nk = _validate_nk(nk, Val{D}())
    _dk = _validate_dk(dk, Val{D}())
    _kmax = _nk .* _dk ./ 2
    return _nk, _kmax
end
function _validate_wavenumber_params(::Nothing, kmax, dk, ::SpatialData{T, D}) where {T, D}
    _kmax = _validate_kmax(kmax, Val{D}())
    _dk = _validate_dk(dk, Val{D}())
    _nk = ceil.(Int, 2 .* _kmax ./ _dk)
    _nk = _validate_nk(_nk, Val{D}())
    return _nk, _kmax
end

function _validate_wavenumber_params(nk::Nothing, kmax::Nothing, dk::Nothing, dim)
    error("Must specify at least two of nk, kmax, or dk.")
end
function _validate_wavenumber_params(nk::Nothing, kmax::Nothing, dk, dim)
    error("Must specify at least two of nk, kmax, or dk, not just dk.")
end
function _validate_wavenumber_params(nk::Nothing, kmax, dk::Nothing, dim)
    error("Must specify at least two of nk, kmax, or dk, not just kmax.")
end
function _validate_wavenumber_params(nk, kmax::Nothing, dk::Nothing, dim)
    error("Must specify at least two of nk, kmax, or dk, not just nk.")
end

# default tapers
function _validate_tapers(tapers::TaperFamily, region, nw)
    return tapers
end
function _validate_tapers(::Nothing, region::Box, nw)
    @argcheck all(nw .> 0)
    M = floor(Int, 2 .* nw .+ 1)
    M = _validate_dims(M, Val{embeddim(region)}())
    return sin_taper_family(M, region)
end
function _validate_tapers(::Nothing, region, nw)
    @argcheck nw isa Int
    @argcheck nw > 0
    bandwidth = nw / Meshes.ustrip(minimum(sides(boundingbox(region))))
    return make_tapers(region, bandwidth = bandwidth)[1]
end
