struct CenteredLFunction{R, L, I, T, D, P, Q} <: IsotropicEstimate{D, P, Q}
    radii::R
    centered_L_function::L
    processinformation::I
    estimationinformation::T
    function CenteredLFunction(
            radii::R, centered_L_function::L, processinfo::ProcessInformation{D},
            estimationinfo::T) where {R, L, T, D}
        P, Q = checkinputs(radii, centered_L_function, processinfo)
        new{R, L, typeof(processinfo), T, D, P, Q}(
            radii, centered_L_function, processinfo, estimationinfo)
    end
end

getargument(f::CenteredLFunction) = f.radii
getestimate(f::CenteredLFunction) = f.centered_L_function

function L2centeredL(radii, k_function)
    [k .- r for (k, r) in zip(k_function, radii)]
end

function centered_L_function(l::LFunction{R, T, D, P}) where {R, T, D, P}
    return CenteredLFunction(l.radii, L2centeredL(l.radii, l.L_function),
        getprocessinformation(l), getestimationinformation(l))
end

centered_L_function(k::KFunction) = centered_L_function(L_function(k))

function centered_L_function(
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
        region,
        radii;
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
        freq_radii = freq_radii,
        rotational_method = rotational_method
    )
    return centered_L_function(k)
end
