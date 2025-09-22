struct PartialKFunction{R, K, I, T, D, P, Q} <: IsotropicEstimate{D, P, Q}
    radii::R
    partial_K_function::K
    processinformation::I
    estimationinformation::T
    function PartialKFunction(
            radii::R, partial_K_function::K, processinfo::ProcessInformation{D},
            estimationinfo::T) where {R, K, T, D}
        P, Q = checkinputs(radii, partial_K_function, processinfo)
        new{R, K, typeof(processinfo), T, D, P, Q}(
            radii, partial_K_function, processinfo, estimationinfo)
    end
end

getargument(f::PartialKFunction) = f.radii
getestimate(f::PartialKFunction) = f.partial_K_function

function partial_K_function(c::PartialCFunction{
        R, C, I, T, D, P, Q}) where {R, C, I, T, D, P, Q}
    mean_prod = getprocessinformation(c).mean_product
    return PartialKFunction(
        c.radii, C2K(c.radii, c.partial_C_function, mean_prod, Val{D}()),
        getprocessinformation(c), getestimationinformation(c))
end

function partial_K_function(spec::PartialSpectra; radii)
    partial_K_function(partial_C_function(spec; radii = radii))
end

function partial_K_function(
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
    c = partial_C_function(
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
    return partial_K_function(c)
end
