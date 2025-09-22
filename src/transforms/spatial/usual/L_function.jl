struct LFunction{R, L, I, T, D, P, Q} <: IsotropicEstimate{D, P, Q}
    radii::R
    L_function::L
    processinformation::I
    estimationinformation::T
    function LFunction(radii::R, L_function::L, processinfo::ProcessInformation{D},
            estimationinfo::T) where {R, L, T, D}
        P, Q = checkinputs(radii, L_function, processinfo)
        new{R, L, typeof(processinfo), T, D, P, Q}(
            radii, L_function, processinfo, estimationinfo)
    end
end

getargument(f::LFunction) = f.radii
getestimate(f::LFunction) = f.L_function

function K2L(k_function, ::Val{D}) where {D}
    V = unitless_measure(Ball(Point(ntuple(x -> 0, Val{D}())), 1))
    return [sign.(k) .* (abs.(k) ./ V) .^ (1 / D) for k in k_function] # broadcast because may contain SVectors
end

function K2L(k_function, ::Val{2})
    return [sign.(k) .* sqrt.(abs.(k) ./ pi) for k in k_function] # broadcast because may contain SVectors
end

function L_function(k::KFunction{R, K, I, T, D, P, Q}) where {R, K, I, T, D, P, Q}
    return LFunction(k.radii, K2L(k.K_function, Val{D}()),
        getprocessinformation(k), getestimationinformation(k))
end

function L_function(
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
    k = K_function(
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
    return L_function(k)
end
