struct CFunction{E, D, P, Q, A, T, IP, IE} <: IsotropicEstimate{E, D, P, Q}
    radii::A
    value::T
    processinformation::IP
    estimationinformation::IE
    function CFunction{E}(radii::A, value::T, processinfo::ProcessInformation{D},
            estimationinfo::IE) where {E, A, T, IE, D}
        P, Q = checkinputs(radii, value, processinfo)
        IP = typeof(processinfo)
        new{E, D, P, Q, A, T, IP, IE}(radii, value, processinfo, estimationinfo)
    end
end
getbaseestimatename(::Type{<:CFunction}) = "C function"
getargument(f::CFunction) = f.radii
getestimate(f::CFunction) = f.value

c_function(data, region; kwargs...) = c_function(spatial_data(data, region); kwargs...)
function c_function(
        data::SpatialData; radii, nfreq, fmax,
        freq_radii = default_rotational_radii(nfreq, fmax),
        rotational_method = default_rotational_kernel(nfreq, fmax), spectra_kwargs...)
    spectrum = spectra(data; nfreq, fmax, spectra_kwargs...)
    return c_function(spectrum; radii = radii, freq_radii = freq_radii,
        rotational_method = rotational_method)
end

function c_function(
        spectrum::Spectra; radii, freq_radii = default_rotational_radii(spectrum),
        rotational_method = default_rotational_kernel(spectrum))
    _c_function(spectrum, radii, freq_radii, rotational_method)
end

function _c_function(spectrum::Spectra{E}, radii, freq_radii, rotational_method) where {E}
    rot_spec = rotational_estimate(spectrum, radii = freq_radii, kernel = rotational_method) # just returns f_mt if rotational_method = NoRotational()
    value = sdf2C(rot_spec, radii)
    return CFunction{E}(
        radii, value, getprocessinformation(spectrum), getestimationinformation(spectrum))
end

function c_function(spectrum::RotationalSpectra{E}; radii) where {E}
    value = sdf2C(rot_spec, radii)
    return CFunction{E}(
        radii, value, getprocessinformation(spectrum), getestimationinformation(spectrum))
end

function partial_c_function(data, region; kwargs...)
    partial_c_function(spatial_data(data, region); kwargs...)
end
function partial_c_function(
        data::SpatialData; radii, nfreq, fmax, spectra_kwargs...)
    f_mt = partial_spectra(data; nfreq, fmax, spectra_kwargs...)
    return c_function(f_mt; radii = radii)
end

function partial_c_function(
        spectrum::NormalOrRotationalSpectra{PartialTrait}; kwargs...)
    c_function(spectrum; kwargs...)
end
function partial_c_function(
        spectrum::NormalOrRotationalSpectra{MarginalTrait};
        kwargs...)
    c_function(partial_spectra(spectrum); kwargs...)
end

## internals

"""
    sdf2C(f, radii::AbstractVector{<:Number})

Takes some form of spectra and returns the C function for each radius in `radii`.
"""
function sdf2C(f, radii::AbstractVector{<:Number})
    [_sdf2C(f, radius) for radius in radii]
end

"""
    _sdf2C(f, radius::Number)

Takes some form of spectra and returns the C function for the `radius`.
"""
function _sdf2C(f::Spectra, radius::Number)
    freq = getargument(f)
    spectra = getestimate(f)
    zero_atom = getprocessinformation(f).atoms
    sdf2C_aniso(freq, spectra, process_trait(f), zero_atom, radius)
end

function sdf2C_aniso(freq, power, zero_atom, radii::AbstractVector{<:Number})
    [sdf2C_aniso(freq, power, zero_atom, radius) for radius in radii]
end
function sdf2C_aniso(
        freq::NTuple{D}, power::AbstractArray{<:Number, D}, ::SingleProcessTrait,
        zero_atom, radii::Number) where {D}
    _sdf2C_aniso(freq, power, zero_atom, radii)
end
function sdf2C_aniso(
        freq::NTuple{D}, power::AbstractArray{<:SMatrix, D}, ::MultipleTupleTrait,
        zero_atom, radii::Number) where {D}
    _sdf2C_aniso(freq, power, zero_atom, radii)
end
function sdf2C_aniso(
        freq::NTuple{D}, power::AbstractArray{<:Number, N}, ::MultipleVectorTrait,
        zero_atom, radii::Number) where {D, N}
    @argcheck length(freq) <= ndims(power)
    out = mapslices(
        z -> _sdf2C_aniso(freq, z, zero_atom, radii), power; dims = (N - D + 1):ndims(y))
    return reshape(out, size(out)[1:(N - D + 1)])
end

function _sdf2C_aniso(freq, power::AbstractArray, zero_atom, radius::Number)
    prod(step, freq) *
    real(sum((s - zero_atom) * sphere_weight(radius, k, Val{length(freq)}())
    for (s, k) in zip(power, Iterators.product(freq...))))
end
function _sdf2C_aniso(
        freq, power::AbstractArray, ::Nothing, radius::Number)
    prod(step, freq) * real(sum(s * sphere_weight(radius, k, Val{length(freq)}())
    for (s, k) in zip(power, Iterators.product(freq...))))
end

function sphere_weight(r, u, ::Val{1})
    x = norm(u)
    return 2r * sinc(2r * x)
end

function sphere_weight(r, u, ::Val{2})
    x = norm(u)
    return (x < 1e-10) ? (pi * r^2) : ((r / x) * besselj1(2π * r * x))
end

function sphere_weight(r, u, ::Val{D}) where {D}
    x = norm(u)
    return (x < 1e-10) ? (unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), r))) :
           (r / x)^(D / 2) * besselj(D / 2, 2π * r * x)
end

##
"""
    _sdf2C(f, radius::Number)

Takes some form of spectra and returns the C function for the `radius`.
"""
function _sdf2C(f::IsotropicEstimate{E, D, P}, radius::Number) where {E, D, P}
    freq = getargument(f)
    spectra = getestimate(f)
    zero_atom = getprocessinformation(f).atoms
    spacing = step(freq)
    # iso weight essentially contains spacing already, because it integrates over an interval
    real(sum((s - zero_atom) * iso_weight(radius, k, spacing, Val{D}())
    for (s, k) in zip(spectra, freq)))
end

function sdf2C_iso(freq, power, zero_atom, radii::AbstractVector{<:Number})
    [sdf2C_iso(freq, power, zero_atom, radius) for radius in radii]
end
function sdf2C_iso(
        freq::NTuple{D}, power::AbstractArray{<:Number, D}, ::SingleProcessTrait,
        zero_atom, radii::Number) where {D}
    _sdf2C_iso(freq, power, zero_atom, radii, Val{D}())
end
function sdf2C_iso(
        freq::NTuple{D}, power::AbstractArray{<:SMatrix, D}, ::MultipleTupleTrait,
        zero_atom, radii::Number) where {D}
    _sdf2C_iso(freq, power, zero_atom, radii, Val{D}())
end
function sdf2C_iso(
        freq::NTuple{D}, power::AbstractArray{<:Number, N}, ::MultipleVectorTrait,
        zero_atom, radii::Number) where {D, N}
    @argcheck length(freq) <= ndims(power)
    out = mapslices(
        z -> _sdf2C_iso(freq, z, zero_atom, radii, Val{D}()), power; dims = (N - D + 1):ndims(y))
    return reshape(out, size(out)[1:(N - D + 1)])
end

function _sdf2C_iso(freq, power, zero_atom, radius::Number, ::Val{D}) where {D}
    real(sum((s - zero_atom) * iso_weight(radius, k, spacing, Val{D}())
    for (s, k) in zip(power, freq)))
end
function _sdf2C_iso(
        freq, power::AbstractArray, ::Nothing, radius::Number, ::Val{D}) where {D}
    real(sum(s * iso_weight(radius, k, spacing, Val{D}())
    for (s, k) in zip(power, freq)))
end

function iso_weight(r, k, s, ::Val{2})
    halfs = s / 2
    besselj0(2π * r * (k - halfs)) - besselj0(2π * r * (k + halfs))
end
