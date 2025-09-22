module SpatialMultitaper

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
include("estimate_types.jl")
include("tapers.jl")
include("mean.jl")
include("tapered_dft.jl")
include("input_checking.jl")
include("spectral_estimate.jl")
include("atoms.jl")

include("transforms/spectral_matrix_transforms.jl")
include("transforms/wavenumber/usual/complex_coherence.jl")
include("transforms/wavenumber/usual/magnitude_coherence.jl")
include("transforms/wavenumber/usual/magnitude_coherence2.jl")
include("transforms/wavenumber/usual/phase.jl")
include("transforms/wavenumber/usual/rotational.jl")

include("transforms/wavenumber/partial/partial_complex_coherence.jl")
include("transforms/wavenumber/partial/partial_magnitude_coherence.jl")
include("transforms/wavenumber/partial/partial_magnitude_coherence2.jl")
include("transforms/wavenumber/partial/partial_phase.jl")
include("transforms/wavenumber/partial/partial_spectra.jl")

include("transforms/spatial/spatial_utils.jl")
include("transforms/spatial/usual/C_function.jl")
include("transforms/spatial/usual/K_function.jl")
include("transforms/spatial/usual/L_function.jl")
include("transforms/spatial/usual/centered_L_function.jl")
include("transforms/spatial/partial/partial_C_function.jl")
include("transforms/spatial/partial/partial_K_function.jl")
include("transforms/spatial/partial/partial_L_function.jl")
include("transforms/spatial/partial/partial_centered_L_function.jl")
include("transforms/spatial/usual/pair_correlation_function.jl")
include("transforms/spatial/partial/partial_pair_correlation_function.jl")

# include("resampling/partial_null_generation.jl")
include("resampling/partial_cross_resampling.jl")
include("resampling/partial_marginal_resampling.jl")
include("resampling/partial_resampling.jl")
include("resampling/resampling.jl")

export multitaper_estimate,
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
       complex_coherence,
       magnitude_coherence,
       magnitude_coherence2,
       phase,
       rotational_estimate,
       partial_spectra,
       partial_complex_coherence,
       partial_magnitude_coherence,
       partial_magnitude_coherence2,
       partial_phase,
       C_function,
       K_function,
       L_function,
       centered_L_function,
       partial_C_function,
       partial_K_function,
       partial_L_function,
       partial_centered_L_function,
       paircorrelation_function,
       paircorrelation_function_direct,
       partial_paircorrelation_function,
       partial_paircorrelation_function_direct,
       shift_resample,
       ToroidalShift,
       partial_shift_resample,
       make_tapers,
       taper_ft,
       check_tapers_for_data

end
