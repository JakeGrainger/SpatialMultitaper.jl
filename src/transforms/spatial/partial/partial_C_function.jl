struct PartialCFunction{R, C, I, T, D, P, Q} <: IsotropicEstimate{D, P, Q}
    radii::R
    partial_C_function::C
    processinformation::I
    estimationinformation::T
    function PartialCFunction(
            radii::R, partial_C_function::C, processinfo::ProcessInformation{D},
            estimationinfo::T) where {R, C, T, D}
        P, Q = checkinputs(radii, partial_C_function, processinfo)
        new{R, C, typeof(processinfo), T, D, P, Q}(
            radii, partial_C_function, processinfo, estimationinfo)
    end
end

getargument(f::PartialCFunction) = f.radii
getestimate(f::PartialCFunction) = f.partial_C_function

function partial_C_function(f::Union{PartialSpectra, IsotropicEstimate}; radii)
    C = sdf2C(f, radii)
    return PartialCFunction(radii, C, getprocessinformation(f), getestimationinformation(f))
end

function partial_C_function(
        f::SpectralEstimate;
        radii,
        rotational_method = NoRotational(),
        freq_radii = default_rotational_radii(f)
)
    partial_C_function(
        rotational_estimate(
            partial_spectra(f),
            radii = freq_radii,
            kernel = rotational_method
        );
        radii = radii
    )
end

function partial_C_function(
        data::Union{NTuple{P, Union{GeoTable, PointSet}}, GeoTable, PointSet},
        region;
        radii,
        nfreq,
        fmax,
        tapers,
        mean_method::MeanEstimationMethod = DefaultMean(),
        freq_radii = default_rotational_radii(nfreq, fmax),
        rotational_method = default_rotational_kernel(freq_radii)
) where {P}
    fhat = multitaper_estimate(
        data,
        region;
        tapers = tapers,
        nfreq = nfreq,
        fmax = fmax,
        mean_method = mean_method
    )
    return partial_C_function(
        fhat;
        radii = radii,
        rotational_method = rotational_method,
        freq_radii = freq_radii
    )
end
