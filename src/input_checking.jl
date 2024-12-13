function check_spatial_data(data)
	_getdims(x) = x isa PointSet ? embeddim(x) : embeddim(domain(x))
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

function check_mean_method(
	mean_method::MeanEstimationMethod,
	data::NTuple{N, Union{GeoTable, PointSet}},
) where {N}
	return ntuple(i -> mean_method, Val{N}())
end

function check_mean_method(
	mean_method::NTuple{P, MeanEstimationMethod},
	data::NTuple{N, Union{GeoTable, PointSet}},
) where {P, N}
	P === N ||
		throw(ArgumentError("Number of mean methods should match number of processes"))
	return mean_method
end

function check_mean_method(
	mean_method::MeanEstimationMethod,
	data::AbstractVector{<:Union{GeoTable, PointSet}},
)
	return fill(mean_method, length(data))
end

function check_mean_method(
	mean_method::AbstractVector{MeanEstimationMethod},
	data::AbstractVector{<:Union{GeoTable, PointSet}},
)
	length(mean_method) == length(data) ||
		throw(ArgumentError("Number of mean methods should match number of processes"))
	return mean_method
end