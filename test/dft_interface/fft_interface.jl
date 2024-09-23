import SpatialMultitaper: pad, choose_freq_res, fft_anydomain, 
    unwrap_fft_output, choose_freq_1d, CartesianGrid, fft, fftshift, make_freq
@testset "fft_interface" begin
    
    @testset "pad" begin
        @test pad(zeros(4,4,4),(10,10,10)) == zeros(10,10,10)
        y = pad(ones(2,2,2),(10,10,10))
        @test all(x->x==1, y[1:2,1:2,1:2])
        @test all(x->x==0, y[3:10,3:10,3:10])
        @test_throws AssertionError pad(ones(3,5), (2,2))
        @test pad(1:3, (5,)) == [1,2,3,0,0]
    end

    @testset "choose_freq_res" begin
        @test choose_freq_res(10,10) == 1        
        @test choose_freq_res(10,12) == 1
        @test choose_freq_res(10,8)  == 2
        @test choose_freq_res(10,4)  == 3
    end
    
    @testset "unwrap_fft_output_1d" begin # nyquist set to n/2 for integer values
        nfreqs = [10,11,30,31,100,101]
        fmaxrels = 1:11
        @testset "nfreq=$nfreq, fmaxrel=$fmaxrel" for nfreq in nfreqs, fmaxrel in fmaxrels
            x = Int.(unwrap_fft_output(choose_freq_1d(nfreq, nfreq/2), (fmaxrel,)))
            y = Int.(choose_freq_1d(nfreq, fmaxrel*(nfreq/2)))
            @test mod.(x,nfreq) ≈ mod.(y,nfreq) atol=1e-4
        end
    end

    @testset "fft_anydomain" begin
        @testset "1d" begin
            @testset "basic" begin
                y = [1,2,3,4]
                y_1_extra = hcat(y, y)
                grid = CartesianGrid((-0.5,), (3.5,), dims=(4,))
                @test fft_anydomain(y, grid, (4,), (1/2,)) ≈ [-2.0, -2.0-2.0im, 10.0, -2.0+2.0im]
                @test fft_anydomain(y, grid, (8,), (1/2,)) ≈ fftshift(fft([1,2,3,4,0,0,0,0])) # 8 is the smallest oversampling
                @test fft_anydomain(y, grid, (8,), (1/2,)) ≈ slow_dft(0:3, y, choose_freq_1d(8, 1/2), -1)

                @test fft_anydomain(y_1_extra, grid, (4,), (1/2,)) ≈ hcat([-2.0, -2.0-2.0im, 10.0, -2.0+2.0im], [-2.0, -2.0-2.0im, 10.0, -2.0+2.0im])

                grid2 = CartesianGrid((1-0.5,), (1+3.5,), dims=(4,))
                @test fft_anydomain(y, grid2, (8,), (1/2,)) ≈ slow_dft(1:4, y, choose_freq_1d(8, 1/2), -1)
            end
            
            starts = [0,2.4]
            y = [3.0, 2.7, 1.8, 2.1, 0.1, -0.5, 0.2, 0.3]
            y_1_extra = hcat(
                [3.0,2.7,1.8,2.1,0.1,-0.5,0.2,0.3],
                [2.0,4.5,1.2,2.3,0.2,-0.1,0.3,0.4],
                [1.0,2.5,1.3,2.5,0.3,-0.2,0.4,0.5],
                [0.0,1.5,1.4,2.7,0.4,-0.3,0.5,0.6]
            )
            y_2_extra = cat(hcat(
                [3.0,2.7,1.8,2.1,0.1,-0.5,0.2,0.3],
                [2.0,4.5,1.2,2.3,0.2,-0.1,0.3,0.4],
                [1.0,2.5,1.3,2.5,0.3,-0.2,0.4,0.5],
                [0.0,1.5,1.4,2.7,0.4,-0.3,0.5,0.6]
            ),hcat(
                [3.1,2.8,1.9,2.0,0.2,-0.6,0.1,0.4],
                [2.1,4.6,1.1,2.4,0.3,-0.2,0.2,0.5],
                [1.1,2.6,1.4,2.6,0.4,-0.1,0.3,0.6],
                [0.1,1.4,1.5,2.8,0.5,-0.4,0.6,0.7]
            ), dims = 3)

            nfreqs = [4,8,10,161]
            fmaxs = 1:8
            @testset "start=$start, nfreq=$nfreq, fmax=$fmax" for start in starts, nfreq in nfreqs, fmax in fmaxs
                grid = CartesianGrid((start-0.25,), (start+4-0.25,), dims=(8,))
                x = start:0.5:start+3.75
                @test fft_anydomain(y, grid, (nfreq,), (fmax,)) ≈ slow_dft(x, y, choose_freq_1d(nfreq, fmax), -1)
                Y_1_extra = fft_anydomain(y_1_extra, grid, (nfreq,), (fmax,))
                for i in axes(Y_1_extra, 2)
                    @test Y_1_extra[:, i] ≈ slow_dft(x, y_1_extra[:, i], choose_freq_1d(nfreq, fmax), -1)
                end
                Y_2_extra = fft_anydomain(y_2_extra, grid, (nfreq,), (fmax,))
                for i in axes(Y_2_extra, 2), j in axes(Y_2_extra, 3)
                    @test Y_2_extra[:, i, j] ≈ slow_dft(x, y_2_extra[:, i, j], choose_freq_1d(nfreq, fmax), -1)
                end
            end
        end

        @testset "2d" begin
            y = [
                3.0 2.7 1.8 2.1 0.1 -0.5 0.2 0.3
                2.0 4.5 1.2 2.3 0.2 -0.1 0.3 0.4
                1.0 2.5 1.3 2.5 0.3 -0.2 0.4 0.5
                0.0 1.5 1.4 2.7 0.4 -0.3 0.5 0.6
            ]
            y_1_extra = cat([
                3.0 2.7 1.8 2.1 0.1 -0.5 0.2 0.3
                2.0 4.5 1.2 2.3 0.2 -0.1 0.3 0.4
                1.0 2.5 1.3 2.5 0.3 -0.2 0.4 0.5
                0.0 1.5 1.4 2.7 0.4 -0.3 0.5 0.6
            ],[
                3.1 2.8 1.9 2.0 0.2 -0.6 0.1 0.4
                2.1 4.6 1.1 2.4 0.3 -0.2 0.2 0.5
                1.1 2.6 1.4 2.6 0.4 -0.1 0.3 0.6
                0.1 1.4 1.5 2.8 0.5 -0.4 0.6 0.7
            ], dims = 3)
            starts = [(0,0), (2.2,2.2), (1.0,-3.4)]
            nfreqs = [(4,4),(8,8),(10,4),(10,15)]
            fmaxs = [(1,1),(1,2), (1,3), (2,2), (2,4), (4,7)]
            @testset "start=$start, nfreq=$nfreq, fmax=$fmax" for start in starts, nfreq in nfreqs, fmax in fmaxs
                grid = CartesianGrid((start[1]-0.25,start[2]-0.25), (start[1]+2-0.25,start[2]+4-0.25), dims=(4,8))
                x = collect(Iterators.product(start[1]:0.5:start[1]+1.75, start[2]:0.5:start[2]+3.75))
                freq = Iterators.ProductIterator(make_freq(nfreq, fmax, 2))
                @test fft_anydomain(y, grid, nfreq, fmax) ≈ slow_dft(x, y, freq, -1)
                Y_1_extra = fft_anydomain(y_1_extra, grid, nfreq, fmax)
                for i in axes(Y_1_extra, 3)
                    @test Y_1_extra[:, :, i] ≈ slow_dft(x, y_1_extra[:, :, i], freq, -1)
                end
            end
        end
    end

end