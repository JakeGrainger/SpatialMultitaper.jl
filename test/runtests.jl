using SpatialMultitaper
using Test

function slow_dft(u, f, freq, iflag)
	pm = iflag ≥ 0 ? 1 : -1
	return [
		sum(f[i] * exp(pm * 2pi * 1im * sum(u[i] .* k)) for i in eachindex(u, f)) for
		k in freq
	]
end
#=
Don't add your tests to runtests.jl. Instead, create files named

    test-title-for-my-test.jl

The file will be automatically included inside a `@testset` with title "Title For My Test".
=#
function get_tests(root)
    for (root, dirs, files) in walkdir(root)
        for file in files
            if isnothing(match(r"^test-.*\.jl$", file))
                continue
            end
            title = titlecase(replace(splitext(file[6:end])[1], "-" => " "))
            @testset "$title" begin
                include(file)
            end
        end
        for dir in dirs
            title = titlecase(replace(dir, "-" => " "))
            @testset begin title
                get_tests(joinpath(root, dir))
            end
        end
    end
end
get_tests(@__DIR__)