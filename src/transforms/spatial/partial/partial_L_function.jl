struct PartialLFunction{R, L, I, T, D, P, Q} <: IsotropicEstimate{D, P, Q}
    radii::R
    partial_L_function::L
    processinformation::I
    estimationinformation::T
    function PartialLFunction(
            radii::R, partial_L_function::L, processinfo::ProcessInformation{D},
            estimationinfo::T) where {R, L, T, D}
        P, Q = checkinputs(radii, partial_L_function, processinfo)
        new{R, L, typeof(processinfo), T, D, P, Q}(
            radii, partial_L_function, processinfo, estimationinfo)
    end
end

getargument(f::PartialLFunction) = f.radii
getestimate(f::PartialLFunction) = f.partial_L_function

function partial_L_function(k::PartialKFunction{
        R, K, I, T, D, P, Q}) where {R, K, I, T, D, P, Q}
    return PartialLFunction(
        k.radii, K2L(k.partial_K_function, Val{D}()), getprocessinformation(k),
        getestimationinformation(k))
end

function partial_L_function(spec::PartialSpectra; radii)
    partial_L_function(partial_K_function(spec; radii = radii))
end

function partial_L_function(
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
    k = partial_K_function(
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
    return partial_L_function(k)
end
