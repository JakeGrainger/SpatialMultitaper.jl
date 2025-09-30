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
    @argcheck all(@. nk isa Integer)
    @argcheck all(@. nk > 0)
    return _validate_dims(nk, Val{D}())
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
