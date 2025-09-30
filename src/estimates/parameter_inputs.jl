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
function _validate_dims(x::NTuple, dim)
    @argcheck length(x) == dim
    return x
end
function _validate_dims(x::Number, dim)
    return ntuple(Returns(x), dim)
end

function _validate_nk(nk, dim)
    @argcheck all(@. nk isa Integer)
    @argcheck all(@. nk > 0)
    return _validate_dims(nk, dim)
end

function _validate_dk(dk, dim)
    @argcheck all(@. dk isa Number)
    @argcheck all(@. dk > 0)
    return _validate_dims(dk, dim)
end

function _validate_kmax(kmax, dim)
    @argcheck all(@. kmax isa Number)
    @argcheck all(@. kmax > 0)
    return _validate_dims(kmax, dim)
end

function _validate_wavenumber_params(nk, kmax, dk::Nothing, dim)
    _nk = _validate_nk(nk, dim)
    _kmax = _validate_kmax(kmax, dim)
    return _nk, _kmax
end
function _validate_wavenumber_params(nk, kmax::Nothing, dk, dim)
    _nk = _validate_nk(nk, dim)
    _dk = _validate_dk(dk, dim)
    _kmax = _nk .* _dk ./ 2
    return _nk, _kmax
end
function _validate_wavenumber_params(nk::Nothing, kmax, dk, dim)
    _kmax = _validate_kmax(kmax, dim)
    _dk = _validate_dk(dk, dim)
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
