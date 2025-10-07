module RCallExt

using SpatialMultitaper
using RCall

import SpatialMultitaper: IsotropicEstimate, SpatialData, getargument, getestimate,
                          getestimatename
import RCall.CategoricalArrays: CategoricalArray, levels
import RCall: rcopy, RClass, rcopytype, sexp, protect, unprotect, setclass!, sexpclass,
              setattrib!

## Isotropic estimates to fv
sexpclass(::IsotropicEstimate) = RClass{:fv}

function sexp(::Type{RClass{:fv}}, est::IsotropicEstimate)
    r = protect(sexp(Dict(:r => getargument(est), :values => getestimate(est))))
    setattrib!(r, "argu", "r")
    setattrib!(r, "valu", "values")
    setattrib!(r, "ylab", getestimatename(est))
    setattrib!(r, "yexp", getestimatename(est))
    setattrib!(r, "fmla", "r~values")
    setattrib!(r, "alim", [extrema(getestimate(est))...])
    setattrib!(r, "labl", ["r", getestimatename(est)])
    setattrib!(r, "desc", ["r", getestimatename(est)])
    setattrib!(r, "units", "NULL")
    setattrib!(r, "fname", "NULL")
    setattrib!(r, "dotnames", "NULL")
    setattrib!(r, "shade", "NULL")

    setclass!(r, sexp("fv"))
    unprotect(1)
    return r
end

## ppp to spatial data

rcopytype(::Type{RClass{:ppp}}, s::Ptr{VecSxp}) = SpatialData

function rcopy(::Type{SpatialData}, s::Ptr{VecSxp})
    spatstat2data(s)
end

function spatstat2data(points)
    data, names = spatstat2points(points)
    region = spatstat2window(points)
    return spatial_data(data, region, names)
end

function spatstat2points(points)
    x = rcopy(points[:x])
    y = rcopy(points[:y])
    markformat = rcopy(points[:markformat])
    if markformat == "none"
        return PointSet(Point.(x, y)), 1
    end
    marks = rcopy(points[:marks])
    if marks isa CategoricalArray
        names = Array(levels(marks))
        points = [PointSet([Point(x[i], y[i]) for i in eachindex(x) if marks[i] == level])
                  for level in names]
        if length(points) <= 10
            return (points...,), names
        else
            return points, names
        end
        return points, names
    elseif marks isa Vector{<:Number}
        return georef((mark = marks,), PointSet(Point.(x, y))), 1
    else
        error("Marks of $(typeof(marks)), format $(markformat) not supported")
    end
end

function spatstat2window(points)
    win_type = rcopy(points[:window][:type])
    if win_type == "rectangle"
        x = rcopy(points[:window][:xrange])
        y = rcopy(points[:window][:yrange])
        return Box(Point(x[1], y[1]), Point(x[2], y[2]))
    elseif win_type == "polygonal"
        bdry = rcopy(points[:window][:bdry])
        bdry_tuples = [[(x, y) for (x, y) in zip(b[:x], b[:y])] for b in bdry]
        if length(bdry_tuples) == 1
            return PolyArea(bdry_tuples[1])
        else
            return PolyArea(bdry_tuples)
        end
    else
        error("Window type $(points[:window][:type]) not supported")
    end
end

end
