struct KFunction{R, K, I, T, D, P, Q} <: IsotropicEstimate{D, P, Q}
    radii::R
    K_function::K
    processinformation::I
    estimationinformation::T
    function KFunction(radii::R, K_function::K, processinfo::ProcessInformation{D},
            estimationinfo::T) where {R, K, T, D}
        P, Q = checkinputs(radii, K_function, processinfo)
        new{R, K, typeof(processinfo), T, D, P, Q}(
            radii, K_function, processinfo, estimationinfo)
    end
end

getargument(f::KFunction) = f.radii
getestimate(f::KFunction) = f.K_function

function C2K(radii, c_function, mean_prod, ::Val{D}) where {D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    return [c ./ mean_prod .+ (V * (r^D)) for (r, c) in zip(radii, c_function)]
end

function K_function(c::CFunction{R, C, I, T, D, P, Q}) where {R, C, I, T, D, P, Q}
    mean_prod = getprocessinformation(c).mean_product
    return KFunction(c.radii, C2K(c.radii, c.C_function, mean_prod, Val{D}()),
        getprocessinformation(c), getestimationinformation(c))
end

function K_function(
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
    c = C_function(
        data,
        region;
        radii = radii,
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
        freq_radii = freq_radii,
        rotational_method = rotational_method
    )
    return K_function(c)
end
