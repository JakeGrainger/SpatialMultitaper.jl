using SpatialMultitaper
using Test

## utilities
function slow_dft(u,f,freq,iflag)
    pm = iflag ≥ 0 ? 1 : -1
    return [sum(f[i]*exp(pm*2pi*1im*sum(u[i].*k)) for i in eachindex(u,f)) for k in freq]
end

@testset "SpatialMultitaper.jl" begin
    include("utils.jl")
    include("tapers.jl")
    
    include("dft_interface/frequencies.jl")
    include("dft_interface/nufft_interface.jl")
    include("dft_interface/fft_interface.jl")

    include("mean.jl")
    include("api.jl")
    include("partial_covariance_density.jl")
    include("general_tapers.jl")

    include("SlepianSolver/runtests.jl")
    include("spectral_matrix_transforms.jl")
    include("K_function.jl")
    include("resampling.jl")
end
