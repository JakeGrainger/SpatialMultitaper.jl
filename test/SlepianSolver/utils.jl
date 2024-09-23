import SpatialMultitaper.upsample
@testset "utils.jl" begin
    @testset "downsample" begin
        x = 1:10
        @test downsample(x, 2) == 1:2:9
        @test downsample(x, 2) !== 2:2:10
        @test_throws AssertionError downsample(x, 11)
        @test_throws AssertionError downsample(x, 0)
        y = rand(10,10,10)
        @test downsample(y, 2) == y[1:2:end,1:2:end,1:2:end]
        @test downsample(y, nothing) === y
    end
    @testset "pad" begin
        M = ones(10,10)
        @test pad(M, 100) == [[ones(10,10); zeros(90,10)] zeros(100,90)]
        @test pad(1:4, 100) == [1:4; zeros(96)]
        @test_throws AssertionError pad(M, (5,11))
    end
    @testset "grid2side" begin
        g = CartesianGrid((100,100), (1.0,2.0), (0.1, 0.3))
        grid = grid2side(g)
        @test length(grid[1]) == 100
        @test length(grid[2]) == 100
        @test grid[1][1] == 1.05
        @test grid[2][1] == 2.15
        @test step.(grid) == (0.1, 0.3)
    end
    @testset "upsample" begin
        @testset "1d" begin
            g = CartesianGrid((100,), (1.0,), (0.1,))
            x = rand(10,)
            @test upsample(x, g, nothing) == x
            @test_throws AssertionError upsample(x,g,2)
            @test_throws AssertionError upsample(x,g,(2,))
            @test size(upsample(x,g,10)) == (100,)
            @test downsample(upsample(x,g,10),10) ≈ x
        end
        @testset "2d" begin
            g = CartesianGrid((100,100), (1.0,2.0), (0.1, 0.3))
            x = rand(10,10)
            @test upsample(x, g, nothing) == x
            @test_throws AssertionError upsample(x,g,2)
            @test_throws AssertionError upsample(x,g,(2,10))
            @test size(upsample(x,g,10)) == (100,100)
            @test downsample(upsample(x,g,10),10) ≈ x
        end
    end
end