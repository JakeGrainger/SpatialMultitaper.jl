module SlepianSolver

using FFTW, LinearAlgebra, Arpack, Meshes, Interpolations
import Base: size, eltype

include("spaces.jl")
include("operators.jl")
include("diagonalisation.jl")

export make_concentration_operator, compute_eigenfunctions

end
