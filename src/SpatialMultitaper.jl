module SpatialMultitaper

using Reexport
@reexport using Meshes, GeoTables

using FFTW,
    FINUFFT,
    Interpolations,
    InvertedIndices,
    LinearAlgebra,
    SpecialFunctions,
    StaticArrays,
    StatsBase

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
include("spectral_estimate.jl")
include("spectral_matrix_transforms.jl")
include("partial_covariance_density.jl")
include("K_function.jl")
include("resampling.jl")

export multitaper_estimate,
    sin_taper_family,
    interpolated_taper_family,
    partial_covariance_density,
    DefaultMean,
    KnownMean,
    georef,
    Point,
    CartesianGrid,
    PointSet,
    pad,
    downsample,
    grid2side,
    side2grid,
    optimaltapers,
    complex_coherence,
    magnitude_coherence,
    magnitude_sq_coherence,
    group_delay,
    partial_spectra,
    partial_complex_coherence,
    partial_magnitude_coherence,
    partial_magnitude_sq_coherence,
    partial_group_delay,
    partial_K,
    shift_resample,
    ToroidalShift,
    partial_K_resample

end
