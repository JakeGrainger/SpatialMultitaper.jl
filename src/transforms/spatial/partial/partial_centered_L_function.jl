struct PartialCenteredLFunction{R, L, I, T, D, P, Q} <: IsotropicEstimate{D, P, Q}
    radii::R
    partial_centered_L_function::L
    processinformation::I
    estimationinformation::T
    function PartialCenteredLFunction(
            radii::R, partial_centered_L_function::L, processinfo::ProcessInformation{D},
            estimationinfo::T) where {R, L, T, D}
        P, Q = checkinputs(radii, partial_centered_L_function, processinfo)
        new{R, L, typeof(processinfo), T, D, P, Q}(
            radii, partial_centered_L_function, processinfo, estimationinfo)
    end
end

getargument(f::PartialCenteredLFunction) = f.radii
getestimate(f::PartialCenteredLFunction) = f.partial_centered_L_function

function partial_centered_L_function(l::PartialLFunction{R, T, D, P}) where {R, T, D, P}
    PartialCenteredLFunction(l.radii, L2centeredL(l.radii, l.L_function),
        getprocessinformation(l), getestimationinformation(l))
end

function partial_centered_L_function(k::PartialKFunction)
    partial_centered_L_function(partial_L_function(k))
end

function partial_centered_L_function(spec::PartialSpectra; radii)
    partial_centered_L_function(partial_K_function(spec; radii = radii))
end

function partial_centered_L_function(
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
    return partial_centered_L_function(k)
end
