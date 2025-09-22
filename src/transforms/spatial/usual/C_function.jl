struct CFunction{R, C, I, T, D, P, Q} <: IsotropicEstimate{D, P, Q}
    radii::R
    C_function::C
    processinformation::I
    estimationinformation::T
    function CFunction(radii::R, C_function::C, processinfo::ProcessInformation{D},
            estimationinfo::T) where {R, C, T, D}
        P, Q = checkinputs(radii, C_function, processinfo)
        new{R, C, typeof(processinfo), T, D, P, Q}(
            radii, C_function, processinfo, estimationinfo)
    end
end

getargument(f::CFunction) = f.radii
getestimate(f::CFunction) = f.C_function

function C_function(f::AbstractEstimate{D, P, N}; radii) where {D, P, N}
    C = sdf2C(f, radii)
    return CFunction(
        radii, C, getprocessinformation(f), getestimationinformation(f))
end

function C_function(
        data,
        region;
        radii,
        nfreq,
        fmax,
        tapers,
        mean_method::MeanEstimationMethod = DefaultMean(),
        freq_radii = default_rotational_radii(nfreq, fmax),
        rotational_method = default_rotational_kernel(freq_radii)
)
    f_mt = multitaper_estimate(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        mean_method = mean_method
    )
    f = rotational_estimate(f_mt, radii = freq_radii, kernel = rotational_method) # just returns f_mt if NoRotational()
    return C_function(f, radii = radii)
end

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
function _sdf2C(
        f::Union{SpectralEstimate{D, F, P, N}, PartialSpectra{D, F, P, N}},
        radius::Number
) where {D, F, P, N}
    freq = getargument(f)
    spectra = getestimate(f)
    zeroatom = getprocessinformation(f).atoms
    _sdf2C(freq, spectra, zeroatom, radius)
end

function _sdf2C(
        freq,
        power::AbstractArray{<:SArray},
        zero_atom,
        radii::AbstractVector{<:Number}
)
    [_sdf2C(freq, power, zero_atom, radius) for radius in radii]
end
function _sdf2C(freq, power::AbstractArray{<:SArray}, zero_atom, radius::Number)
    prod(step, freq) * real(
        sum(
        (s - zero_atom) * sphere_weight(radius, k, Val{length(freq)}())
    for
    (s, k) in zip(power, Iterators.product(freq...))
    ),
    )
end

function _sdf2C(freq, power::AbstractArray{<:SArray}, zero_atom::Nothing, radius::Number)
    prod(step, freq) * real(
        sum(
        s * sphere_weight(radius, k, Val{length(freq)}())
    for
    (s, k) in zip(power, Iterators.product(freq...))
    ),
    )
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
function _sdf2C(f::IsotropicEstimate{D, P}, radius::Number) where {D, P}
    freq = getargument(f)
    spectra = getestimate(f)
    zero_atom = getprocessinformation(f).atoms
    spacing = step(freq)
    real(
        sum(
        (s - zero_atom) * iso_weight(radius, k, spacing, Val{D}())
    for
    (s, k) in zip(spectra, freq)
    ),
    )
end

function iso_weight(r, k, s, ::Val{2})
    halfs = s / 2
    besselj0(2π * r * (k - halfs)) - besselj0(2π * r * (k + halfs))
end
