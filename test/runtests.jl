using SpatialMultitaper
import SpatialMultitaper as Spmt
using Test
include("SpatialMultitaperTestingUtils.jl")
using .SpatialMultitaperTestingUtils

#=
Don't add your tests to runtests.jl. Instead, create files named

    test-title-for-my-test.jl

The file will be automatically included inside a `@testset` with title "Title For My Test".
=#
function get_tests(folder)
    for (root, dirs, files) in walkdir(folder)
        for dir in dirs
            title = titlecase(replace(dir, "-" => " "))
            @testset "$title" begin
                get_tests(joinpath(root, dir))
            end
        end
        for file in files
            if isnothing(match(r"^test-.*\.jl$", file))
                continue
            end
            title = titlecase(replace(splitext(file[6:end])[1], "-" => " "))
            @testset "$title" begin
                include(joinpath(root, file))
            end
        end
    end
end
get_tests(@__DIR__)