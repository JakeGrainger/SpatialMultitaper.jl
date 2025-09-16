using SpatialMultitaper, Test, SafeTestsets

@safetestset "utils" begin
    include("test-utils.jl")
end

#=
Don't add your tests to runtests.jl. Instead, create files named

    test-title-for-my-test.jl

The file will be automatically included inside a `@testset` with title "Title For My Test".
=#

include("SpatialMultitaperTestingUtils.jl")
using .SpatialMultitaperTestingUtils
import SpatialMultitaper as Spmt
for (root, dirs, files) in walkdir(@__DIR__)
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
