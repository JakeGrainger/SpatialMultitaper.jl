using SpatialMultitaper, Test
import SpatialMultitaper: _choose_wavenumbers_1d

@testset "wavenumbers" begin
    @test _choose_wavenumbers_1d(5, 2.5) == -2.0:2.0
    @test _choose_wavenumbers_1d(6, 3) == -3.0:2.0
    @test _choose_wavenumbers_1d(5, 1 / 2) == (-2.0:2.0) ./ 5
    @test _choose_wavenumbers_1d(6, 1 / 2) == (-3.0:2.0) ./ 6
end
