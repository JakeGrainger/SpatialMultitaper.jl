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
    ConstructionBase

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

include("transforms/spectral_matrix_transforms.jl")
include("transforms/usual/complex_coherence.jl")
include("transforms/usual/magnitude_coherence.jl")
include("transforms/usual/magnitude_coherence2.jl")
include("transforms/usual/phase.jl")
include("transforms/partial/partial_complex_coherence.jl")
include("transforms/partial/partial_magnitude_coherence.jl")
include("transforms/partial/partial_magnitude_coherence2.jl")
include("transforms/partial/partial_phase.jl")
include("transforms/partial/partial_spectra.jl")

include("K_function.jl")
include("resampling.jl")

export multitaper_estimate,
    sin_taper_family,
    interpolated_taper_family,
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
    magnitude_coherence2,
    phase,
    partial_spectra,
    partial_complex_coherence,
    partial_magnitude_coherence,
    partial_magnitude_coherence2,
    partial_phase,
    partial_K,
    shift_resample,
    ToroidalShift,
    partial_K_resample,
    make_tapers,
    taper_ft,
    check_tapers_for_data

end
