using Test, SpatialMultitaper, StableRNGs
include("../../test_utilities/TestData.jl")
using .TestData

import SpatialMultitaper: getestimate

# Create test data
rng = StableRNG(123)
data = make_points_example(rng, n_processes = 3, dim = 2,
    return_type = :vector, point_number = 100)

@testset "rotational_spectra" begin
    result = rotational_spectra(data, nk = 100, kmax = 1.0)
    @test all(isreal.(getestimate(result)))

    # With custom radii
    radii = range(0, 0.5, length = 20)
    result_custom = rotational_spectra(data, nk = 100, kmax = 1.0, radii = radii)
    @test size(getestimate(result_custom), 3) == length(radii)
end

@testset "rotational_partial_spectra" begin
    result = rotational_partial_spectra(data, nk = 100, kmax = 1.0)
    @test all(isreal.(getestimate(result)))
end

@testset "rotational_coherence" begin
    result = rotational_coherence(data, nk = 100, kmax = 1.0)
    @test all(isreal.(getestimate(result)))
    @test all(-1 .<= getestimate(result) .<= 1)
end

@testset "rotational_partial_coherence" begin
    result = rotational_partial_coherence(data, nk = 100, kmax = 1.0)
    @test all(isreal.(getestimate(result)))
    @test all(-1 .<= getestimate(result) .<= 1)
end
