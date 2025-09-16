using SpatialMultitaper, Test

import SpatialMultitaper: unitless_coords, unitless_spacing, unitless_origin,
                          unitless_minimum, unitless_maximum, unitless_measure,
                          point2unitlesstype, points2coords, box2sides, grid2sides,
                          padto, centerpad, downsample, upsample, mask

"""
    basic_test_grid(T, min, max, dims)

Create a test `CartesianGrid` with specified type `T`, min, max, and dims.
"""
function basic_test_grid(T, min, max, dims)
    return CartesianGrid(convert.(T, min), convert.(T, max), dims = dims)
end

@testset "unitless_coords" begin
    types = (Float32, Float64)
    @testset "testing types $T" for T in types
        point1d = Point(convert(T, 1))
        point2d = Point(convert(T, 1), convert(T, 2))
        point3d = Point(convert(T, 1), convert(T, 2), convert(T, 3))

        @test unitless_coords(point1d) == (1,)
        @test unitless_coords(point2d) == (1, 2)
        @test unitless_coords(point3d) == (1, 2, 3)
        @test unitless_coords(point1d) isa Tuple{T}
        @test unitless_coords(point2d) isa Tuple{T, T}
        @test unitless_coords(point3d) isa Tuple{T, T, T}
    end
end

@testset "unitless_spacing" begin
    types = (Float32, Float64)
    @testset "testing types $T" for T in types
        grid1d = basic_test_grid(T, (1,), (2,), (10,))
        grid2d = basic_test_grid(T, (1, 2), (3, 4), (10, 10))
        grid3d = basic_test_grid(T, (1, 2, 3), (4, 5, 6), (10, 10, 10))

        @test unitless_spacing(grid1d) == convert.(T, ((2 - 1) / 10,))
        @test unitless_spacing(grid2d) == convert.(T, ((3 - 1) / 10, (4 - 2) / 10))
        @test unitless_spacing(grid3d) ==
              convert.(T, ((4 - 1) / 10, (5 - 2) / 10, (6 - 3) / 10))
        @test unitless_spacing(grid1d) isa Tuple{T}
        @test unitless_spacing(grid2d) isa Tuple{T, T}
        @test unitless_spacing(grid3d) isa Tuple{T, T, T}
    end
end

@testset "unitless_origin" begin
    types = (Float32, Float64)
    @testset "testing types $T" for T in types
        grid1d = basic_test_grid(T, (1,), (2,), (10,))
        grid2d = basic_test_grid(T, (1, 2), (3, 4), (10, 10))
        grid3d = basic_test_grid(T, (1, 2, 3), (4, 5, 6), (10, 10, 10))

        @test unitless_origin(grid1d) == convert.(T, (1,))
        @test unitless_origin(grid2d) == convert.(T, (1, 2))
        @test unitless_origin(grid3d) == convert.(T, (1, 2, 3))
        @test unitless_origin(grid1d) isa Tuple{T}
        @test unitless_origin(grid2d) isa Tuple{T, T}
        @test unitless_origin(grid3d) isa Tuple{T, T, T}
    end
end

@testset "unitless_minimum" begin
    types = (Float32, Float64)
    @testset "testing types $T" for T in types
        grid1d = basic_test_grid(T, (1,), (2,), (10,))
        grid2d = basic_test_grid(T, (1, 2), (3, 4), (10, 10))
        grid3d = basic_test_grid(T, (1, 2, 3), (4, 5, 6), (10, 10, 10))

        @test unitless_minimum(grid1d) == convert.(T, (1,))
        @test unitless_minimum(grid2d) == convert.(T, (1, 2))
        @test unitless_minimum(grid3d) == convert.(T, (1, 2, 3))
        @test unitless_minimum(grid1d) isa Tuple{T}
        @test unitless_minimum(grid2d) isa Tuple{T, T}
        @test unitless_minimum(grid3d) isa Tuple{T, T, T}
    end
end

@testset "unitless_maximum" begin
    types = (Float32, Float64)
    @testset "testing types $T" for T in types
        grid1d = basic_test_grid(T, (1,), (2,), (10,))
        grid2d = basic_test_grid(T, (1, 2), (3, 4), (10, 10))
        grid3d = basic_test_grid(T, (1, 2, 3), (4, 5, 6), (10, 10, 10))

        @test unitless_maximum(grid1d) == convert.(T, (2,))
        @test unitless_maximum(grid2d) == convert.(T, (3, 4))
        @test unitless_maximum(grid3d) == convert.(T, (4, 5, 6))
        @test unitless_maximum(grid1d) isa Tuple{T}
        @test unitless_maximum(grid2d) isa Tuple{T, T}
        @test unitless_maximum(grid3d) isa Tuple{T, T, T}
    end
end

@testset "unitless_measure" begin
    types = (Float32, Float64)
    @testset "testing types $T" for T in types
        box1d = Box(Point(convert(T, 1)), Point(convert(T, 2)))
        box2d = Box(Point(convert(T, 1), convert(T, 2)),
            Point(convert(T, 3), convert(T, 4)))
        box3d = Box(
            Point(convert(T, 1), convert(T, 2), convert(T, 3)),
            Point(convert(T, 4), convert(T, 5), convert(T, 6))
        )

        @test unitless_measure(box1d) == convert(T, 1)
        @test unitless_measure(box2d) == convert(T, 4)
        @test unitless_measure(box3d) == convert(T, 27)
        @test unitless_measure(box1d) isa T
        @test unitless_measure(box2d) isa T
        @test unitless_measure(box3d) isa T
    end
end

@testset "point2unitlesstype" begin
    types = (Float32, Float64)
    @testset "testing types $T" for T in types
        point1d = Point(convert(T, 1))
        point2d = Point(convert(T, 1), convert(T, 2))
        point3d = Point(convert(T, 1), convert(T, 2), convert(T, 3))

        @test point2unitlesstype(point1d) == T
        @test point2unitlesstype(point2d) == T
        @test point2unitlesstype(point3d) == T
    end
end

@testset "points2coords" begin
    types = (Float32, Float64)
    @testset "testing types $T" for T in types
        points1d = PointSet([Point(convert(T, i)) for i in 1:10])
        points2d = PointSet([Point(convert(T, i), convert(T, i + 1)) for i in 1:10])
        points3d = PointSet([Point(convert(T, i), convert(T, i + 1),
                                 convert(T, i + 2)) for i in 1:10])

        coords1d = points2coords(points1d)
        coords2d = points2coords(points2d)
        coords3d = points2coords(points3d)

        @test coords1d == (convert.(T, 1:10),)
        @test coords2d == (convert.(T, 1:10), convert.(T, 2:11))
        @test coords3d == (convert.(T, 1:10), convert.(T, 2:11), convert.(T, 3:12))

        @test coords1d[1] isa AbstractVector{T}
        @test coords2d[1] isa AbstractVector{T}
        @test coords2d[2] isa AbstractVector{T}
        @test coords3d[1] isa AbstractVector{T}
        @test coords3d[2] isa AbstractVector{T}
        @test coords3d[3] isa AbstractVector{T}
    end
end

@testset "box2sides" begin
    types = (Float32, Float64)
    @testset "testing types $T" for T in types
        box1d = Box(Point(convert(T, 1)), Point(convert(T, 2)))
        box2d = Box(Point(convert(T, 1), convert(T, 2)),
            Point(convert(T, 3), convert(T, 4)))
        box3d = Box(
            Point(convert(T, 1), convert(T, 2), convert(T, 3)),
            Point(convert(T, 4), convert(T, 5), convert(T, 6))
        )
        @test box2sides(box1d) == (convert.(T, (1, 2)),)
        @test box2sides(box2d) ==
              (convert.(T, (1, 3)), convert.(T, (2, 4)))
        @test box2sides(box3d) == (
            convert.(T, (1, 4)),
            convert.(T, (2, 5)),
            convert.(T, (3, 6))
        )
        @test box2sides(box1d) isa NTuple{1, Tuple{T, T}}
        @test box2sides(box2d) isa NTuple{2, Tuple{T, T}}
        @test box2sides(box3d) isa NTuple{3, Tuple{T, T}}
    end
end

@testset "grid2sides" begin
    types = (Float64,) # TODO: Float32 version currently broken in Meshes.jl
    @testset "testing types $T" for T in types
        grid1d = basic_test_grid(T, (1,), (2,), (10,))
        grid2d = basic_test_grid(T, (1, 2), (3, 4), (10, 10))
        grid3d = basic_test_grid(T, (1, 2, 3), (4, 5, 6), (10, 10, 10))
        gridsides1d = grid2sides(grid1d)
        gridsides2d = grid2sides(grid2d)
        gridsides3d = grid2sides(grid3d)

        @test gridsides1d isa NTuple{1, AbstractVector{T}}
        @test gridsides2d isa NTuple{2, AbstractVector{T}}
        @test gridsides3d isa NTuple{3, AbstractVector{T}}

        @test size(grid1d) == length.(gridsides1d)
        @test size(grid2d) == length.(gridsides2d)
        @test size(grid3d) == length.(gridsides3d)

        @test unitless_spacing(grid1d) == step.(gridsides1d)
        @test unitless_spacing(grid2d) == step.(gridsides2d)
        @test unitless_spacing(grid3d) == step.(gridsides3d)

        @test centroid(grid1d, 1) ≈ Point(getindex.(gridsides1d, 1))
        @test centroid(grid2d, 1) ≈ Point(getindex.(gridsides2d, 1))
        @test centroid(grid3d, 1) ≈ Point(getindex.(gridsides3d, 1))
        @test centroid(grid1d, 10) ≈ Point(getindex.(gridsides1d, 10))
        @test centroid(grid2d, 100) ≈ Point(getindex.(gridsides2d, 10))
        @test centroid(grid3d, 1000) ≈ Point(getindex.(gridsides3d, 10))
    end
end

@testset "padto" begin
    vec = ones(10)
    mat = ones(10, 10)
    array = ones(10, 10, 10)
    @test padto(vec, 15) == [ones(10); zeros(5)]
    @test padto(mat, (15, 20)) == [ones(10, 10) zeros(10, 10); zeros(5, 20)]
    array_padded = padto(array, (15, 20, 25))
    @test all(
        idx -> checkbounds(Bool, array, idx) ? array_padded[idx] == 1.0 :
               array_padded[idx] == 0.0, CartesianIndices(array_padded))

    @test_throws ArgumentError padto(mat, (2, 2))
    @test_throws ArgumentError padto(mat, (100, 2))
    @test_throws MethodError padto(vec, (2, 2))

    @test eltype(padto(1:4, 10)) == Int
    @test eltype(padto(1.0f0:4.0f0, 10)) == Float32
    @test eltype(padto(1.0:4.0, 10)) == Float64
end

@testset "centerpad" begin
    vec = 1:4
    mat = reshape(1:9, 3, 3)
    array = reshape(1:27, 3, 3, 3)
    @test centerpad(vec, 1) == [0; 1:4; 0]
    @test centerpad(mat, 1) == [0 0 0 0 0
                                0 1 4 7 0
                                0 2 5 8 0
                                0 3 6 9 0
                                0 0 0 0 0]
    array_padded = centerpad(array, 1)
    @test array_padded[2:4, 2:4, 2:4] == array
    @test_throws ArgumentError centerpad(vec, -1)
    @test_throws ArgumentError centerpad(mat, (1, -1))
    @test_throws ArgumentError centerpad(vec, 1.0)
end

@testset "downsample" begin
    vec = 1:20
    mat = reshape(1:100, 10, 10)
    array = reshape(1:1000, 10, 10, 10)
    @test downsample(vec, 2) == vec[1:2:end]
    @test downsample(mat, (2, 3)) == mat[1:2:end, 1:3:end]
    @test downsample(array, (2, 3, 4)) == array[1:2:end, 1:3:end, 1:4:end]
    @test downsample(vec, nothing) == vec
    @test downsample(mat, nothing) == mat
    @test downsample(array, nothing) == array

    @test_throws ArgumentError downsample(vec, 0)
    @test_throws ArgumentError downsample(mat, (2, 0))
    @test_throws ArgumentError downsample(mat, (2, -1))
    @test_throws ArgumentError downsample(mat, (11, 2))
    @test_throws ArgumentError downsample(mat, (2, 11))
    @test_throws MethodError downsample(vec, (2, 3))

    @test eltype(downsample(1:20, 2)) == Int
    @test eltype(downsample(1.0f0:20.0f0, 2)) == Float32
    @test eltype(downsample(1.0:20.0, 2)) == Float64
end

@testset "upsample" begin
    @testset "1d" begin
        g = CartesianGrid((100,), (1.0,), (0.1,))
        x = rand(10)
        @test upsample(x, g, nothing) == x
        @test size(upsample(x, g, 10)) == (100,)
        @test downsample(upsample(x, g, 10), 10) ≈ x

        @test_throws DimensionMismatch upsample(x, g, 2)
        @test_throws DimensionMismatch upsample(x, g, (2,))
        @test_throws DimensionMismatch upsample(randn(10, 10), g, 2)
    end
    @testset "2d" begin
        g = CartesianGrid((100, 100), (1.0, 2.0), (0.1, 0.3))
        x = rand(10, 10)
        @test upsample(x, g, nothing) == x
        @test size(upsample(x, g, 10)) == (100, 100)
        @test downsample(upsample(x, g, 10), 10) ≈ x

        @test_throws DimensionMismatch upsample(x, g, 2)
        @test_throws DimensionMismatch upsample(x, g, (2, 10))
    end
end

@testset "mask" begin
    @testset "points" begin
        points = PointSet(vec([Point(i, j) for i in 1:5, j in 1:5]))
        region = Box(Point(1.5, 1.5), Point(4.5, 4.5))
        masked_points = mask(points, region)
        @test masked_points == PointSet(vec([Point(i, j) for i in 2:4, j in 2:4]))
        @test all(p ∈ region for p in masked_points)

        marked_points = georef((mark = 1:25,), points)
        masked_marked_points = mask(marked_points, region)
        surviving_marks = reshape(marked_points.mark, 5, 5)[2:4, 2:4][:]
        @test masked_marked_points == georef(
            (mark = surviving_marks,), PointSet(vec([Point(i, j) for i in 2:4, j in 2:4])))
    end

    @testset "geotables" begin
        grid = CartesianGrid((0, 0), (5, 5), dims = (5, 2))
        geotable = georef((rf = randn(10),), grid)
        region = Box(Point(1.5, 1.5), Point(4.5, 4.5))
        masked_geotable = mask(geotable, region)
        @test all(centroid(grid, i) ∈ region && !isnan(m) || isnan(m)
        for (i, m) in enumerate(masked_geotable.rf))
    end
end
