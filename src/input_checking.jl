_getdims(x) = x isa PointSet ? embeddim(x) : embeddim(domain(x))
function check_spatial_data(data)
    dim = _getdims(first(data))
    @assert all(_getdims.(data) .== dim) "data should all be the same spatial dimension"
    copied_data = deepcopy(data)
    for proc in copied_data
        if proc isa GeoTable
            if length(values(proc)) > 1
                @warn "more than one random field provided to a geotable, currently we only process the first of these!"
            end
            if any(x -> abs(x) == Inf, values(proc)[1])
                error("Some fields have infinite values!")
            end
            if any(isnan, values(proc)[1])
                replace!(values(proc)[1], NaN => zero(eltype(values(proc)[1])))
            end
        end
    end
    return copied_data, dim
end
