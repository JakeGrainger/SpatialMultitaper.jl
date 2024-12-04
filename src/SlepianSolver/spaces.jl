abstract type Space end
struct TensorSpace{T, D} <: Space
	m::NTuple{D, Int}
end
Base.size(s::TensorSpace) = s.m
Base.size(s::TensorSpace{T, 1}) where {T} = s.m[1]
Base.eltype(::TensorSpace{T, D}) where {T, D} = T

struct AnySpace <: Space end
Base.size(::AnySpace) = spaceerror()
Base.eltype(::AnySpace) = spaceerror()
spaceerror() = error("This operator is not specific enough, compose with identity to fix.")

struct AnyElementSpace{D} <: Space
	m::NTuple{D, Int}
end
Base.size(s::AnyElementSpace) = s.m
Base.size(s::AnyElementSpace{1}) = s.m[1]
Base.eltype(::AnyElementSpace) = spaceerror()

struct AnySizeSpace{T} <: Space end
Base.size(::AnySizeSpace) = spaceerror()
Base.eltype(::AnySizeSpace{T}) where {T} = T

function check_size(t1::Space, t2::Space)
	@assert size(t1) == size(t2)
	nothing
end

check_size(::Union{AnySpace, AnySizeSpace}, ::Union{AnySpace, AnySizeSpace}) = nothing
check_size(::Union{AnySpace, AnySizeSpace}, ::Space) = nothing
check_size(::Space, ::Union{AnySpace, AnySizeSpace}) = nothing

function check_eltype(t1::Space, t2::Space)
	@assert eltype(t1) == eltype(t2)
	nothing
end

check_eltype(::Union{AnySpace, AnyElementSpace}, ::Union{AnySpace, AnyElementSpace}) =
	nothing
check_eltype(::Union{AnySpace, AnyElementSpace}, ::Space) = nothing
check_eltype(::Space, ::Union{AnySpace, AnyElementSpace}) = nothing
