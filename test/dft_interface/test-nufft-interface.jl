using SpatialMultitaper, Test
import SpatialMultitaper: _choose_wavenumbers_1d, nufft1d1_anydomain, nufft2d1_anydomain,
                          nufft3d1_anydomain, rescale_points, freq_downsample_startindex,
                          box2sides, points2coords, nufft_anydomain
include("../test_utilities/TestUtils.jl")
using .TestUtils: slow_dft

function test_nufft1d1_anydomain(interval, nk, kmax, xj, cj)
    freq = _choose_wavenumbers_1d(nk, kmax)
    fast_out = nufft1d1_anydomain(interval, nk, kmax, xj, cj, -1, 1e-14)[:, 1]
    slow_out = slow_dft(xj, cj, freq, -1)
    fast_out_pos = nufft1d1_anydomain(interval, nk, kmax, xj, cj, 1, 1e-14)[:, 1]
    slow_out_pos = slow_dft(xj, cj, freq, 1)
    @test length(fast_out) == length(slow_out)
    @test fast_out ≈ slow_out
    @test fast_out_pos ≈ slow_out_pos
end

function test_nufft2d1_anydomain(box, nk, kmax, xj, yj, cj)
    freq_x = _choose_wavenumbers_1d(nk[1], kmax[1])
    freq_y = _choose_wavenumbers_1d(nk[2], kmax[2])
    freq = Iterators.product(freq_x, freq_y)
    uj = [(xj[i], yj[i]) for i in eachindex(xj, yj)]

    fast_out = nufft2d1_anydomain(box, nk, kmax, xj, yj, cj, -1, 1e-14)[:, :, 1]
    slow_out = slow_dft(uj, cj, freq, -1)
    fast_out_pos = nufft2d1_anydomain(box, nk, kmax, xj, yj, cj, 1, 1e-14)[:, :, 1]
    slow_out_pos = slow_dft(uj, cj, freq, 1)
    @test length(fast_out) == length(slow_out)
    @test fast_out ≈ slow_out
    @test fast_out_pos ≈ slow_out_pos
end

function test_nufft3d1_anydomain(box, nk, kmax, xj, yj, zj, cj)
    freq_x = _choose_wavenumbers_1d(nk[1], kmax[1])
    freq_y = _choose_wavenumbers_1d(nk[2], kmax[2])
    freq_z = _choose_wavenumbers_1d(nk[3], kmax[3])
    freq = Iterators.product(freq_x, freq_y, freq_z)
    uj = [(xj[i], yj[i], zj[i]) for i in eachindex(xj, yj, zj)]

    fast_out = nufft3d1_anydomain(box, nk, kmax, xj, yj, zj, cj, -1, 1e-14)[:, :, :, 1]
    slow_out = slow_dft(uj, cj, freq, -1)
    fast_out_pos = nufft3d1_anydomain(box, nk, kmax, xj, yj, zj, cj, 1, 1e-14)[
        :, :, :, 1]
    slow_out_pos = slow_dft(uj, cj, freq, 1)
    @test length(fast_out) == length(slow_out)
    @test fast_out ≈ slow_out
    @test fast_out_pos ≈ slow_out_pos
end

@testset "rescale_points" begin
    interval = [-1 / 2, 1 / 2]
    nks = [10, 9]
    kmaxs = [5, 4.5]
    x = [0.1, 0.2, 0.3, 0.4, 0.5]
    for (nk, kmax) in zip(nks, kmaxs)
        x_scaled, oversample, shift = rescale_points(x, nk, kmax, interval)
        @test x_scaled == [0.1, 0.2, 0.3, 0.4, 0.5] .* 2π
        @test oversample == 1
        @test rescale_points(x, nk, kmax - 1, interval)[2] == 1
        @test rescale_points(x, nk, kmax + 1, interval)[2] == 2
    end
end

@testset "_choose_wavenumbers_1d" begin
    nks = [9, 11, 20]
    kmax = 0.5
    oversample = [1, 2, 3, 7, 21, 22]
    @testset "nk=$nk, ovesample=$c" for nk in nks, c in oversample
        @test _choose_wavenumbers_1d(nk * c, kmax)[freq_downsample_startindex(
            nk,
            c
        ):c:end] ≈ _choose_wavenumbers_1d(nk, kmax)
    end
end

@testset "1d" begin
    # simple case
    intervals = [[-1 / 2, 1 / 2], [-2.3, 2.3], [-2.3, 1.4]]
    nks = [3, 10, 9, 21]
    kmaxs = [4, 5, 6, 0.6, 1.1]
    xj = [0.24, 0.13, 0.41, 0.43, 0.01]
    cj = complex.([0.4, 0.5, 0.6, 0.7, 0.8])
    @testset "kmax=$kmax, interval=$interval, nk=$nk" for kmax in kmaxs,
        interval in intervals,
        nk in nks

        test_nufft1d1_anydomain(interval, nk, kmax, xj, cj)
    end
end

@testset "2d" begin
    # simple case
    intervals = [
        ([-1 / 2, 1 / 2], [-1 / 2, 1 / 2]),
        ([-2.3, 2.3], [-1.2, 1.2]),
        ([-2.3, 1.4], [-0.6, 4.1])
    ]
    nks = [(3, 4), (10, 8), (9, 11), (21, 22)]
    kmaxs = [(4, 4), (5, 5), (6.4, 0.6), (0.6, 0.7)]
    xj = [0.24, 0.13, 0.41, 0.43, 0.01]
    yj = [0.36, 0.23, 0.51, 0.53, 0.11]
    cj = complex.([0.4, 0.5, 0.6, 0.7, 0.8])
    @testset "kmax=$kmax, interval=$interval, nk=$nk" for kmax in kmaxs,
        interval in intervals,
        nk in nks

        test_nufft2d1_anydomain(interval, nk, kmax, xj, yj, cj)
    end
end

@testset "3d" begin
    # simple case
    intervals = [
        ([-1 / 2, 1 / 2], [-1 / 2, 1 / 2], [-1 / 2, 1 / 2]),
        ([-2.3, 2.3], [-1.2, 1.2], [-6.4, 6.4]),
        ([-2.3, 1.4], [-0.6, 4.1], [-0.8, 1.1])
    ]
    nks = [(3, 4, 5), (10, 8, 9), (9, 11, 12), (21, 22, 23)]
    kmaxs = [(4, 4, 4), (5, 5, 5), (6.4, 0.6, 0.7), (0.6, 0.7, 0.8)]
    xj = [0.24, 0.13, 0.41, 0.43, 0.01]
    yj = [0.36, 0.23, 0.51, 0.53, 0.11]
    zj = [0.76, 0.33, 0.21, 0.76, 0.01]
    cj = complex.([0.4, 0.5, 0.6, 0.7, 0.8])
    @testset "kmax=$kmax, interval=$interval, nk=$nk" for kmax in kmaxs,
        interval in intervals,
        nk in nks

        test_nufft3d1_anydomain(interval, nk, kmax, xj, yj, zj, cj)
    end
end

@testset "box2sides" begin
    @test box2sides(Box(Point(4), Point(7))) == ((4, 7),)
    @test box2sides(Box(Point(0, 3), Point(4, 5))) == ((0, 4), (3, 5))
    @test box2sides(Box(Point(0, 3, -1), Point(4, 5, -0.5))) ==
          ((0, 4), (3, 5), (-1, -0.5))
end

@testset "points2coords" begin
    @test points2coords(PointSet([Point(1, 2), Point(3, 5)])) == ([1, 3], [2, 5])
end

@testset "general" begin
    points1 = PointSet([Point(1), Point(3)])
    points2 = PointSet([Point(1, 2), Point(3, 5)])
    points3 = PointSet([Point(1, 2, 3), Point(3, 5, 7)])
    cj = complex.([3.4, 2.5])
    nk = (4, 4, 4)
    kmax = (1, 1, 1)
    nufft_anydomain(
        Box(Point(0), Point(2)),
        nk[1:1],
        kmax[1:1],
        points1,
        cj,
        -1,
        1e-9
    ) == nufft1d1_anydomain((0, 2), nk[1], kmax[1], [1, 3], cj, -1, 1e-9)
    nufft_anydomain(
        Box(Point(0, 1), Point(2, 5)),
        nk[1:2],
        kmax[1:2],
        points2,
        cj,
        -1,
        1e-9
    ) == nufft2d1_anydomain(
        ((0, 2), (1, 5)),
        nk[1:2],
        kmax[1:2],
        [1, 3],
        [3, 5],
        cj,
        -1,
        1e-9
    )
    nufft_anydomain(
        Box(Point(0, 1, 2), Point(2, 5, 7)),
        nk,
        kmax,
        points3,
        cj,
        -1,
        1e-9
    ) == nufft3d1_anydomain(
        ((0, 2), (1, 5), (2, 7)),
        nk,
        kmax,
        [1, 3],
        [3, 5],
        [5, 7],
        cj,
        -1,
        1e-9
    )
end
