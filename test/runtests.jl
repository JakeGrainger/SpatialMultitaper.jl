using SpatialMultitaper, Test, SafeTestsets, Aqua

# Aqua.test_all(SpatialMultitaper)
run_all = true
run_spatial = true

# Core functionality tests
if run_all
    @safetestset "Covariance zero atom" begin
        include("test-covariance_zero_atom.jl")
    end

    @safetestset "Utils" begin
        include("test-utils.jl")
    end
    @safetestset "General Tapers" begin
        include("tapers/test-general-tapers.jl")
    end
    @safetestset "Tapers" begin
        include("tapers/test-tapers.jl")
    end
    @safetestset "Input Checking" begin
        include("test-input-data.jl")
    end
    @safetestset "Mean" begin
        include("test-mean.jl")
    end

    # # DFT Interface tests
    @safetestset "FFT Interface" begin
        include("dft_interface/test-fft-interface.jl")
    end
    @safetestset "Wavenumbers" begin
        include("dft_interface/test-wavenumbers.jl")
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

    # @safetestset "Parameter Inputs" begin
    #     include("estimates/test-parameter-inputs.jl")
    # end

    @safetestset "Spectral Estimate" begin
        include("estimates/wavenumber/test-spectra.jl")
    end

    @safetestset "Transforms" begin
        include("estimates/wavenumber/test-spectral-matrix-transforms.jl")
    end

    @safetestset "Partial spectra" begin
        include("estimates/wavenumber/test-partial-spectra.jl")
    end

    @safetestset "coherence" begin
        include("estimates/wavenumber/test-coherence.jl")
    end

    @safetestset "Rotational functions" begin
        include("estimates/wavenumber/test-rotational-functions.jl")
    end

    @safetestset "Rotational estimate" begin
        include("estimates/generic/test-rotational-estimate.jl")
    end

    @safetestset "Marginal transform" begin
        include("estimates/generic/test-marginal-transform.jl")
    end

    @safetestset "Printing" begin
        include("estimates/test-printing.jl")
    end

    # @safetestset "Resampling" begin
    #     include("test-resampling.jl")
    # end
    @warn "Skipping resampling tests temporarily"
end
if run_spatial || run_all
    # spatial
    @safetestset "c function" begin
        include("estimates/spatial/test-c-function.jl")
    end

    @safetestset "k function" begin
        include("estimates/spatial/test-k-function.jl")
    end

    @safetestset "L function" begin
        include("estimates/spatial/test-l-function.jl")
    end

    @safetestset "Centered L function" begin
        include("estimates/spatial/test-centered-l-function.jl")
    end
end
