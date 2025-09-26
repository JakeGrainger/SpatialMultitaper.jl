module SpatialMultitaper

# TODO: should use ArgCheck.jl

using Reexport
@reexport using Meshes, GeoTables

using Distributions,
      FFTW,
      FINUFFT,
      Interpolations,
      InvertedIndices,
      LinearAlgebra,
      SpecialFunctions,
      StaticArrays,
      StatsBase,
      ConstructionBase,
      LazyArrays,
      HypergeometricFunctions,
      GeoStatsProcesses,
      Random

import BSplineKit

include("SlepianSolver/SlepianSolver.jl")
using .SlepianSolver

include("general_tapers.jl")

include("dft_interface/frequencies.jl")
include("dft_interface/nufft_interface.jl")
include("dft_interface/fft_interface.jl")

include("utils.jl")
include("tapers.jl")
include("mean.jl")
include("tapered_dft.jl")
include("input_checking.jl")
include("covariance_zero_atom.jl")

include("estimates/estimate_types.jl")
include("estimates/printing.jl")
include("estimates/errors.jl")
include("estimates/rotational_estimate.jl")
include("estimates/marginal_transform.jl")

include("estimates/spectral_estimate.jl")
include("estimates/partial_spectra.jl")

include("estimates/spectral_matrix_transforms.jl")
include("estimates/coherence.jl")

include("estimates/spatial/c_function.jl")
include("estimates/spatial/k_function.jl")
include("estimates/spatial/l_function.jl")
include("estimates/spatial/centered_l_function.jl")
include("estimates/spatial/pair_correlation_function.jl")

# include("resampling/partial_null_generation.jl")
include("resampling/partial_cross_resampling.jl")
include("resampling/partial_marginal_resampling.jl")
include("resampling/partial_resampling.jl")
include("resampling/resampling.jl")

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
       paircorrelation_function_direct,
       partial_paircorrelation_function,
       partial_paircorrelation_function_direct,
       shift_resample,
       ToroidalShift,
       partial_shift_resample,
       make_tapers,
       taper_ft,
       check_tapers_for_data,
       mask

end
