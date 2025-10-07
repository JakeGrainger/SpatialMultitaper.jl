function _spatstat2points(points)
    x = points[:x]
    y = points[:y]
    if points[:markformat] == "none"
        return PointSet(Point.(x, y))
    elseif points[:marks] isa CategoricalArray
        m = points[:marks]
        return ([PointSet([Point(x[i], y[i]) for i in eachindex(x) if m[i] == level])
                 for level in levels(m)]...,)
    else
        error("Marks of $(typeof(points[:marks])), format $(points[:markformat]) not supported")
    end
end
