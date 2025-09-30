using SpatialMultitaper, Test
include("../test_utilities/TestUtils.jl")
using .TestUtils

import SpatialMultitaper: choose_freq_oversample, unwrap_fft_output, _choose_frequencies_1d,
                          fft_anydomain,
                          _make_frequency_grid, unwrap_index, fftshift, fft

@testset "choose_freq_oversample" begin
    # n = 10, nk = 10
    @test choose_freq_oversample(10, 10) == 1
    # n = 10, nk = 12
    @test choose_freq_oversample(10, 12) == 1
    # n = 10, nk = 8
    @test choose_freq_oversample(10, 8) == 2
    # n = 10, nk = 4
    @test choose_freq_oversample(10, 4) == 3
    # nk = 0 (should throw)
    @test_throws ArgumentError choose_freq_oversample(10, 0)
    # n = 0 (should throw)
    @test_throws ArgumentError choose_freq_oversample(0, 10)
    # n = 1, nk = 1
    @test choose_freq_oversample(1, 1) == 1
    # n = 1, nk = 2
    @test choose_freq_oversample(1, 2) == 1
    # Large n, small nk
    @test choose_freq_oversample(1000, 1; maxoversample = 1000) == 1000
    # Large n, large nk
    @test choose_freq_oversample(1000, 1000) == 1
    # nk > n
    @test choose_freq_oversample(10, 20) == 1
    # maxitr reached (should throw)
    @test_throws ErrorException choose_freq_oversample(100, 1; maxoversample = 50)
end

@testset "unwrap_fft_output_1d" begin # nyquist set to n/2 for integer values
    @testset "1d" begin
        # Empty array
        @test_throws ArgumentError unwrap_fft_output(Float64[], (1,))
        # Error: mismatched kmaxrel
        @test_throws ArgumentError unwrap_fft_output(ones(10), (1, 2))
        # Error: negative kmaxrel
        @test_throws ArgumentError unwrap_fft_output(ones(10), (-1,))
        # correctness on frequencies
        nks = [10, 11, 30, 31, 100, 101]
        kmaxrels = 1:11
        @testset "nk=$nk, kmaxrel=$kmaxrel" for nk in nks, kmaxrel in kmaxrels
            # obtained from applying unwrap_fft_output to the fft frequencies spaced by 1
            x = Int.(unwrap_fft_output(
                _choose_frequencies_1d(nk, nk / 2), (kmaxrel,)))
            # the fft frequencies spaced by kmaxrel (which should be the output up to periodicity)
            y = Int.(_choose_frequencies_1d(nk, kmaxrel * (nk / 2)))
            @test mod.(x, nk)≈mod.(y, nk) atol=1e-4
        end
    end

    @testset "2d" begin
        # Error: mismatched kmaxrel
        @test_throws ArgumentError unwrap_fft_output(ones(10, 10), (1, 1, 1))
        # Error: negative kmaxrel
        @test_throws ArgumentError unwrap_fft_output(ones(10, 10), (1, -1))
    end

    @testset "3d" begin
        # Error: mismatched kmaxrel
        @test_throws ArgumentError unwrap_fft_output(ones(10, 10, 10), (1, 1, 1, 1))
        # Error: negative kmaxrel
        @test_throws ArgumentError unwrap_fft_output(ones(10, 10, 10), (1, 1, -1))
    end
end

@testset "unwrap_index" begin
    # Even n, even a
    @test unwrap_index(10, 2) == 6:2:24
    # Even n, odd a
    @test unwrap_index(10, 3) == 1:3:30
    # Odd n, even a
    @test unwrap_index(11, 2) == 7:2:27
    # Odd n, odd a
    @test unwrap_index(11, 3) == 2:3:33
    # a == 1
    @test unwrap_index(10, 1) == 1:1:10
    @test unwrap_index(11, 1) == 1:1:11
    # Edge case: n = 1
    @test unwrap_index(1, 2) == 2:2:2
    # Edge case: a = 0 (should error)
    @test_throws ArgumentError unwrap_index(10, 0)
end

@testset "fft_anydomain" begin
    @testset "1d" begin
        @testset "basic" begin
            y = [1, 2, 3, 4]
            y_1_extra = hcat(y, y)
            grid = CartesianGrid((-0.5,), (3.5,), dims = (4,))
            @test fft_anydomain(y, grid, (4,), (1 / 2,)) ≈
                  [-2.0, -2.0 - 2.0im, 10.0, -2.0 + 2.0im]
            @test fft_anydomain(y, grid, (8,), (1 / 2,)) ≈
                  fftshift(fft([1, 2, 3, 4, 0, 0, 0, 0])) # 8 is the smallest oversampling
            @test fft_anydomain(y, grid, (8,), (1 / 2,)) ≈
                  slow_dft(0:3, y, _choose_frequencies_1d(8, 1 / 2), -1)

            @test fft_anydomain(y_1_extra, grid, (4,), (1 / 2,)) ≈ hcat(
                [-2.0, -2.0 - 2.0im, 10.0, -2.0 + 2.0im],
                [-2.0, -2.0 - 2.0im, 10.0, -2.0 + 2.0im]
            )

            grid2 = CartesianGrid((1 - 0.5,), (1 + 3.5,), dims = (4,))
            @test fft_anydomain(y, grid2, (8,), (1 / 2,)) ≈
                  slow_dft(1:4, y, _choose_frequencies_1d(8, 1 / 2), -1)
        end

        starts = [0, 2.4]
        y = [3.0, 2.7, 1.8, 2.1, 0.1, -0.5, 0.2, 0.3]
        y_1_extra = hcat(
            [3.0, 2.7, 1.8, 2.1, 0.1, -0.5, 0.2, 0.3],
            [2.0, 4.5, 1.2, 2.3, 0.2, -0.1, 0.3, 0.4],
            [1.0, 2.5, 1.3, 2.5, 0.3, -0.2, 0.4, 0.5],
            [0.0, 1.5, 1.4, 2.7, 0.4, -0.3, 0.5, 0.6]
        )
        y_2_extra = cat(
            hcat(
                [3.0, 2.7, 1.8, 2.1, 0.1, -0.5, 0.2, 0.3],
                [2.0, 4.5, 1.2, 2.3, 0.2, -0.1, 0.3, 0.4],
                [1.0, 2.5, 1.3, 2.5, 0.3, -0.2, 0.4, 0.5],
                [0.0, 1.5, 1.4, 2.7, 0.4, -0.3, 0.5, 0.6]
            ),
            hcat(
                [3.1, 2.8, 1.9, 2.0, 0.2, -0.6, 0.1, 0.4],
                [2.1, 4.6, 1.1, 2.4, 0.3, -0.2, 0.2, 0.5],
                [1.1, 2.6, 1.4, 2.6, 0.4, -0.1, 0.3, 0.6],
                [0.1, 1.4, 1.5, 2.8, 0.5, -0.4, 0.6, 0.7]
            ),
            dims = 3
        )

        nks = [4, 8, 10, 161]
        kmaxs = 1:8
        @testset "start=$start, nk=$nk, kmax=$kmax" for start in starts,
            nk in nks,
            kmax in kmaxs

            grid = CartesianGrid((start - 0.25,), (start + 4 - 0.25,), dims = (8,))
            x = start:0.5:(start + 3.75)
            @test fft_anydomain(y, grid, (nk,), (kmax,)) ≈
                  slow_dft(x, y, _choose_frequencies_1d(nk, kmax), -1)
            Y_1_extra = fft_anydomain(y_1_extra, grid, (nk,), (kmax,))
            for i in axes(Y_1_extra, 2)
                @test Y_1_extra[:, i] ≈
                      slow_dft(x, y_1_extra[:, i], _choose_frequencies_1d(nk, kmax), -1)
            end
            Y_2_extra = fft_anydomain(y_2_extra, grid, (nk,), (kmax,))
            for i in axes(Y_2_extra, 2), j in axes(Y_2_extra, 3)
                @test Y_2_extra[:, i, j] ≈
                      slow_dft(
                    x, y_2_extra[:, i, j], _choose_frequencies_1d(nk, kmax), -1)
            end
        end
    end

    @testset "2d" begin
        y = [3.0 2.7 1.8 2.1 0.1 -0.5 0.2 0.3
             2.0 4.5 1.2 2.3 0.2 -0.1 0.3 0.4
             1.0 2.5 1.3 2.5 0.3 -0.2 0.4 0.5
             0.0 1.5 1.4 2.7 0.4 -0.3 0.5 0.6]
        y_1_extra = cat(
            [3.0 2.7 1.8 2.1 0.1 -0.5 0.2 0.3
             2.0 4.5 1.2 2.3 0.2 -0.1 0.3 0.4
             1.0 2.5 1.3 2.5 0.3 -0.2 0.4 0.5
             0.0 1.5 1.4 2.7 0.4 -0.3 0.5 0.6],
            [3.1 2.8 1.9 2.0 0.2 -0.6 0.1 0.4
             2.1 4.6 1.1 2.4 0.3 -0.2 0.2 0.5
             1.1 2.6 1.4 2.6 0.4 -0.1 0.3 0.6
             0.1 1.4 1.5 2.8 0.5 -0.4 0.6 0.7],
            dims = 3
        )
        starts = [(0, 0), (2.2, 2.2), (1.0, -3.4)]
        nks = [(4, 4), (8, 8), (10, 4), (10, 15)]
        kmaxs = [(1, 1), (1, 2), (1, 3), (2, 2), (2, 4), (4, 7)]
        @testset "start=$start, nk=$nk, kmax=$kmax" for start in starts,
            nk in nks,
            kmax in kmaxs

            grid = CartesianGrid(
                (start[1] - 0.25, start[2] - 0.25),
                (start[1] + 2 - 0.25, start[2] + 4 - 0.25),
                dims = (4, 8)
            )
            x = collect(
                Iterators.product(
                start[1]:0.5:(start[1] + 1.75), start[2]:0.5:(start[2] + 3.75)),
            )
            freq = Iterators.ProductIterator(_make_frequency_grid(nk, kmax, 2))
            @test fft_anydomain(y, grid, nk, kmax) ≈ slow_dft(x, y, freq, -1)
            Y_1_extra = fft_anydomain(y_1_extra, grid, nk, kmax)
            for i in axes(Y_1_extra, 3)
                @test Y_1_extra[:, :, i] ≈ slow_dft(x, y_1_extra[:, :, i], freq, -1)
            end
        end
    end
end
