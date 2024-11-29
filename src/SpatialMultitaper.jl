module SpatialMultitaper

using Meshes, FFTW, FINUFFT, Interpolations, StatsBase, GeoTables, LinearAlgebra, SpecialFunctions

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

export multitaper_estimate, sin_taper_family, interpolated_taper_family, partial_covariance_density, DefaultMean, KnownMean
export georef, Point, CartesianGrid, PointSet
export pad, downsample, grid2side, side2grid
export optimaltapers
export complex_coherence, magnitude_coherence, magnitude_sq_coherence, group_delay, 
        partial_complex_coherence, partial_magnitude_coherence, partial_magnitude_sq_coherence, partial_group_delay
export partial_K
export marginal_shift, ToroidalShift

end
