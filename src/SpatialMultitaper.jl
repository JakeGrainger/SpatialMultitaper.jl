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

include("functional_statistics/infrastructure/process_traits.jl")
include("functional_statistics/infrastructure/processinfo.jl")
include("functional_statistics/infrastructure/interface.jl")
include("functional_statistics/infrastructure/default_estimate.jl")
include("functional_statistics/infrastructure/indexing.jl")
include("functional_statistics/infrastructure/names.jl")
include("functional_statistics/infrastructure/misc.jl")
include("functional_statistics/infrastructure/printing.jl")
include("functional_statistics/infrastructure/errors.jl")
include("functional_statistics/infrastructure/spectral_matrix_transforms.jl")

include("functional_statistics/rotational_estimate/rotational_estimate.jl")
include("functional_statistics/marginal_transform/marginal_transform.jl")

include("functional_statistics/spectra/spectra.jl")
include("functional_statistics/partial_spectra/partial_spectra.jl")
include("functional_statistics/coherence/coherence.jl")
include("functional_statistics/rotational_spectra/rotational_spectra.jl")

include("functional_statistics/c_function/c_function.jl")
include("functional_statistics/k_function/k_function.jl")
include("functional_statistics/l_function/l_function.jl")
include("functional_statistics/centered_l_function/centered_l_function.jl")
include("functional_statistics/pair_correlation_function/pair_correlation_function.jl")

include("api/wavenumber/api.jl")
include("api/spatial/api.jl")
include("api/general.jl")

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
