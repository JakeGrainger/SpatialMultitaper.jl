using SpatialMultitaper, Test, SafeTestsets, Aqua

# Aqua.test_all(SpatialMultitaper)

# Core functionality tests
@safetestset "Utils" begin
    include("test-utils.jl")
end
@safetestset "General Tapers" begin
    include("test-general-tapers.jl")
end
@safetestset "Tapers" begin
    include("test-tapers.jl")
end
@safetestset "Input Checking" begin
    include("test-input-checking.jl")
end
@safetestset "Mean" begin
    include("test-mean.jl")
end

# # DFT Interface tests
@safetestset "FFT Interface" begin
    include("dft_interface/test-fft-interface.jl")
end
@safetestset "Frequencies" begin
    include("dft_interface/test-frequencies.jl")
end
@safetestset "NUFFT Interface" begin
    include("dft_interface/test-nufft-interface.jl")
end

# # Slepian Solver tests
@safetestset "Diagonalisation" begin
    include("SlepianSolver/test-diagonalisation.jl")
end
@safetestset "Operator" begin
    include("SlepianSolver/test-operator.jl")
end

# Dft tests
@safetestset "Tapered DFT" begin
    include("test-tapered-dft.jl")
end

# estimate tests
@safetestset "Estimate Types" begin
    include("estimates/test-estimate-types.jl")
end

@safetestset "Error Messages" begin
    include("estimates/test-errors.jl")
end

@safetestset "Spectral Estimate" begin
    include("estimates/test-spectral-estimate.jl")
end

@safetestset "Transforms" begin
    include("estimates/test-spectral-matrix-transforms.jl")
end

@safetestset "Partial spectra" begin
    include("estimates/test-partial-spectra.jl")
end

@safetestset "coherence" begin
    include("estimates/test-coherence.jl")
end

@safetestset "Rotational estimate" begin
    include("estimates/test-rotational-estimate.jl")
end

@safetestset "Marginal transform" begin
    include("estimates/test-marginal-transform.jl")
end

@safetestset "spatial functions" begin
    include("estimates/test-spatial-functions.jl")
end

@safetestset "Printing" begin
    include("estimates/test-printing.jl")
end

@safetestset "Resampling" begin
    include("test-resampling.jl")
end
