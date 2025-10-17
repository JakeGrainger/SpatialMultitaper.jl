module SpatialMultitaper

using Reexport
@reexport using Meshes, GeoTables

using Distributions, FFTW, FINUFFT, Interpolations, InvertedIndices, LinearAlgebra,
      SpecialFunctions, StaticArrays, StatsBase, ConstructionBase, LazyArrays,
      HypergeometricFunctions, GeoStatsProcesses, Random, ArgCheck, CircularArrays

import DataAPI: ncol

import BSplineKit

include("SlepianSolver/SlepianSolver.jl")
using .SlepianSolver

include("utils.jl")

include("dft_interface/wavenumbers.jl")
include("dft_interface/nufft_interface.jl")
include("dft_interface/fft_interface.jl")

include("input_data.jl")
include("tapers/tapers.jl")
include("scalar_statistics/mean.jl")
include("scalar_statistics/covariance_zero_atom.jl")

include("tapered_dft.jl")

include("functional_statistics/process_traits.jl")
include("functional_statistics/processinfo.jl")
include("functional_statistics/estimate_types.jl")
include("functional_statistics/printing.jl")
include("functional_statistics/errors.jl")
include("functional_statistics/parameter_inputs.jl")

include("functional_statistics/generic/rotational_estimate.jl")
include("functional_statistics/generic/marginal_transform.jl")

include("functional_statistics/wavenumber/spectra.jl")
include("functional_statistics/wavenumber/partial_spectra.jl")
include("functional_statistics/wavenumber/spectral_matrix_transforms.jl")
include("functional_statistics/wavenumber/coherence.jl")
include("functional_statistics/wavenumber/rotational_functions.jl")

include("functional_statistics/spatial/c_function.jl")
include("functional_statistics/spatial/k_function.jl")
include("functional_statistics/spatial/l_function.jl")
include("functional_statistics/spatial/centered_l_function.jl")
include("functional_statistics/spatial/pair_correlation_function.jl")

include("resampling/partial_resampling/partial_cross_resampling.jl")
include("resampling/partial_resampling/partial_marginal_resampling.jl")
include("resampling/partial_resampling/partial_resampling.jl")
include("resampling/shift_resampling/shift.jl")

export ncol, observations, getregion, spatial_data

export spectra,
       multitaper_estimate,
       sin_taper_family,
       interpolated_taper_family,
       DefaultMean,
       KnownMean,
       georef,
       Point,
       CartesianGrid,
       PointSet,
       padto,
       downsample,
       optimaltapers,
       coherence,
       magnitude_coherence,
       magnitude_squared_coherence,
       phase,
       rotational_estimate,
       rotational_spectra,
       rotational_partial_spectra,
       rotational_coherence,
       rotational_partial_coherence,
       partial_spectra,
       partial_coherence,
       partial_magnitude_coherence,
       partial_magnitude_squared_coherence,
       partial_phase,
       c_function,
       k_function,
       l_function,
       centered_l_function,
       partial_c_function,
       partial_k_function,
       partial_l_function,
       partial_centered_l_function,
       paircorrelation_function,
       pair_correlation_function_direct,
       partial_paircorrelation_function,
       partial_pair_correlation_function_direct,
       shift_resample,
       ToroidalShift,
       partial_shift_resample,
       make_tapers,
       taper_ft,
       check_tapers_for_data,
       mask

end
