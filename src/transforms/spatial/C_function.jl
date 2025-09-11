struct CFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    C_function::T
    function CFunction(radii::R, C::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, C)
        new{R,T,D,P}(radii, C)
    end
end

getargument(f::CFunction) = f.radii
getestimate(f::CFunction) = f.C_function
getextrafields(::CFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function C_function(f::AbstractEstimate{D,1,N}, zero_atom; radii) where {D,N}
    C = sdf2C(f, zero_atom[1], radii) # one dimensional case, indexing `zero_atom` returns value if a number
    return CFunction(radii, C, Val{D}())
end

function C_function(f::AbstractEstimate{D,P,N}, zero_atom; radii) where {D,P,N}
    C = sdf2C(f, zero_atom, radii)
    return CFunction(radii, C, Val{D}())
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
    rotational_method = default_rotational_kernel(freq_radii),
)
    f_mt = multitaper_estimate(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        mean_method = mean_method,
    )
    f = rotational_estimate(f_mt, radii = freq_radii, kernel = rotational_method) # just returns f_mt if NoRotational()
    zero_atom = atom_estimate(data, region)
    return C_function(f, zero_atom, radii = radii)
end


function sdf2C(f, zero_atom::Number, radii)
    _sdf2C(f, zero_atom, radii)
end

function sdf2C(f, zero_atom::Tuple, radii)
    _sdf2C(f, diagm(SVector(zero_atom...)), radii)
end

"""
    _sdf2C(f, zero_atom, radii::AbstractVector{<:Number})

Takes some form of spectra and returns the C function for each radius in `radii`.
"""
function _sdf2C(f, zero_atom, radii::AbstractVector{<:Number})
    [_sdf2C(f, zero_atom, radius) for radius in radii]
end

"""
    _sdf2C(f, zero_atom, radius::Number)

Takes some form of spectra and returns the C function for the `radius`.
"""
function _sdf2C(
    f::Union{SpectralEstimate{D,F,P,N},PartialSpectra{D,F,P,N}},
    zero_atom,
    radius::Number,
) where {D,F,P,N}
    freq = getargument(f)
    spectra = getestimate(f)
    _sdf2C(freq, spectra, zero_atom, radius)
end

function _sdf2C(
    freq,
    power::AbstractArray{<:SArray},
    zero_atom,
    radii::AbstractVector{<:Number},
)
    [_sdf2C(freq, power, zero_atom, radius) for radius in radii]
end
function _sdf2C(freq, power::AbstractArray{<:SArray}, zero_atom, radius::Number)
    prod(step, freq) * real(
        sum(
            (s - zero_atom) * sphere_weight(radius, k, Val{length(freq)}()) for
            (s, k) in zip(power, Iterators.product(freq...))
        ),
    )
end

function _sdf2C(freq, power::AbstractArray{<:SArray}, zero_atom::Nothing, radius::Number)
    prod(step, freq) * real(
        sum(
            s * sphere_weight(radius, k, Val{length(freq)}()) for
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
    _sdf2C(f, zero_atom, radius::Number)

Takes some form of spectra and returns the C function for the `radius`.
"""
function _sdf2C(f::IsotropicEstimate{D,P}, zero_atom, radius::Number) where {D,P}
    freq = getargument(f)
    spectra = getestimate(f)
    spacing = step(freq)
    real(
        sum(
            (s - zero_atom) * iso_weight(radius, k, spacing, Val{D}()) for
            (s, k) in zip(spectra, freq)
        ),
    )
end

function iso_weight(r, k, s, ::Val{2})
    halfs = s / 2
    besselj0(2π * r * (k - halfs)) - besselj0(2π * r * (k + halfs))
end