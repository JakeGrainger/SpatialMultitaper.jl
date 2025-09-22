using SpatialMultitaper, Test, SafeTestsets

@safetestset "utils" begin
    include("test-utils.jl")
end
@safetestset "fft_interface" begin
    include("dft_interface/test-fft-interface.jl")
end

include("SpatialMultitaperTestingUtils.jl")
using .SpatialMultitaperTestingUtils
import SpatialMultitaper as Spmt
@testset "General Tapers" begin
    include("test-general-tapers.jl")
end
@testset "Input Checking" begin
    include("test-input-checking.jl")
end
@testset "Mean" begin
    include("test-mean.jl")
end
@testset "Resampling" begin
    include("test-resampling.jl")
end
@testset "Spectral Estimate" begin
    include("test-spectral-estimate.jl")
end
@testset "Tapered Dft" begin
    include("test-tapered-dft.jl")
end
@testset "Tapers" begin
    include("test-tapers.jl")
end
@testset "Diagonalisation" begin
    include("SlepianSolver/test-diagonalisation.jl")
end
@testset "Operator" begin
    include("SlepianSolver/test-operator.jl")
end
@testset "Frequencies" begin
    include("dft_interface/test-frequencies.jl")
end
@testset "Nufft Interface" begin
    include("dft_interface/test-nufft-interface.jl")
end
@testset "Partial Complex Coherence" begin
    include("transforms/partial/test-partial-complex-coherence.jl")
end
@testset "Partial Magnitude Coherence" begin
    include("transforms/partial/test-partial-magnitude-coherence.jl")
end
@testset "Partial Magnitude Coherence2" begin
    include("transforms/partial/test-partial-magnitude-coherence2.jl")
end
@testset "Partial Phase" begin
    include("transforms/partial/test-partial-phase.jl")
end
@testset "Partial Spectra" begin
    include("transforms/partial/test-partial-spectra.jl")
end
@testset "C Function" begin
    include("transforms/spatial/test-C-function.jl")
end
@testset "Partial K Function" begin
    include("transforms/spatial/test-partial-K-function.jl")
end
@testset "Complex Coherence" begin
    include("transforms/usual/test-complex-coherence.jl")
end
@testset "Magnitude Coherence" begin
    include("transforms/usual/test-magnitude-coherence.jl")
end
@testset "Magnitude Coherence2" begin
    include("transforms/usual/test-magnitude-coherence2.jl")
end
@testset "Phase" begin
    include("transforms/usual/test-phase.jl")
end
@testset "Rotational" begin
    include("transforms/usual/test-rotational.jl")
end
