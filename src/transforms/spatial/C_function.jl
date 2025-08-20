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

function C_function(
    f::SpectralEstimate{D,F,1,N},
    zero_atom,
    radii,
) where {D,F,N}
    C = sdf2C(f, zero_atom[1], radii) # one dimensional case, indexing `zero_atom` returns value if a number
    return CFunction(radii, C, Val{D}())
end

function C_function(
    f::SpectralEstimate{D,F,P,N},
    zero_atom,
    radii,
) where {D,F,P,N}
    C = sdf2C(f, zero_atom, radii)
    return CFunction(radii, C, Val{D}())
end

function C_function(
    data,
    region,
    radii;
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
)
    f = multitaper_estimate(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        mean_method = mean_method,
    )
    zero_atom = atom_estimate(data, region)
    return C_function(f, zero_atom, radii)
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
    prod(step, freq) * real(
        sum(
            (s - zero_atom) * sphere_weight(radius, k, Val{D}()) for
            (s, k) in zip(spectra, Iterators.product(freq...))
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