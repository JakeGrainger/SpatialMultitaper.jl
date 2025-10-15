# Core type definitions
include("types.jl")

# Taper evaluation and function interface
include("evaluation.jl")

# Construction utilities
include("construction/discrete.jl")
include("construction/interpolated.jl")
include("construction/sin_tapers.jl")
include("construction/general_tapers.jl")

# Validation utilities
include("validation/normalization.jl")
include("validation/concentration.jl")
include("validation/orthogonality.jl")

# Grid interface utilities
include("grid_interface.jl")

# I/O and display methods
include("io.jl")
