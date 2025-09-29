"""
	LinearOperator

A linear operator. Assumes that inputspace and outputspace are defined.
One should also define mul!(y, L::LinearOperator, x) and preallocate_output(L::LinearOperator, x) for efficiency.
"""
abstract type LinearOperator end
Base.eltype(L::LinearOperator) = promote_type(eltype(inputspace(L)), eltype(outputspace(L)))
Base.size(L::LinearOperator) = (size(outputspace(L)), size(inputspace(L)))
Base.size(L::LinearOperator, i::Int) = Base.size(L)[i]
LinearAlgebra.issymmetric(::LinearOperator) = false
Base.:*(L::LinearOperator, x::AbstractArray) = mul!(preallocate_output(L, x), L, x)
preallocate_input(::LinearOperator) = nothing
preallocate_output(::LinearOperator, ::Nothing) = nothing

struct CompositeLinearOperator{L1 <: LinearOperator, L2 <: LinearOperator, A} <:
       LinearOperator
    input::L1
    output::L2
    intermediate::A
end

inputspace(c::CompositeLinearOperator) = inputspace(c.input)
outputspace(c::CompositeLinearOperator) = outputspace(c.output)
function preallocate_output(c::CompositeLinearOperator, x::Nothing) # needed for dispatch
    preallocate_output(c.output, c.intermediate)
end
function preallocate_output(c::CompositeLinearOperator, x)
    preallocate_output(c.output, c.intermediate)
end

function preallocate_output(
        ::CompositeLinearOperator{LinearOperator, LinearOperator, Nothing}, x::Nothing)
    error("Not well specified enough to run composite operator. Please rebuild.")
end
function preallocate_output( # needed for dispatch
        op::CompositeLinearOperator{LinearOperator, LinearOperator, Nothing}, x)
    preallocate_output(op, nothing)
end
preallocate_input(c::CompositeLinearOperator) = preallocate_input(c.input)
function LinearAlgebra.mul!(y::AbstractArray, c::CompositeLinearOperator, x::AbstractArray)
    mul!(c.intermediate, c.input, x)
    mul!(y, c.output, c.intermediate)
    return y
end

function Base.:*(output::LinearOperator, input::LinearOperator)
    check_compatible(input, output)
    _CompositeLinearOperator(
        input,
        output,
        preallocate_output(input, preallocate_input(input))
    )
end

function _CompositeLinearOperator(
        input::LinearOperator, output::LinearOperator, intermediate)
    CompositeLinearOperator(input, output, intermediate)
end
function _CompositeLinearOperator(
        input::LinearOperator,
        output::CompositeLinearOperator{L1, L2, Nothing},
        intermediate
) where {L1, L2}
    if intermediate === nothing
        CompositeLinearOperator(input, output, intermediate)
    else # We rebuild the output intermediate storage if we now know its type from our new input. This will recursively rebuild any missing types.
        rebuilt_output = _CompositeLinearOperator(
            output.input,
            output.output,
            preallocate_output(output.input, intermediate)
        ) # because intermediate is the input type of output.input
        CompositeLinearOperator(input, rebuilt_output, intermediate)
    end
end

function check_compatible(input::LinearOperator, output::LinearOperator)
    check_size(outputspace(input), inputspace(output))
    check_eltype(outputspace(input), inputspace(output))
    nothing
end

struct SupportProjection{T, D, A <: AbstractArray{T, D}} <: LinearOperator
    region::A
end
function LinearAlgebra.mul!(y, s::SupportProjection, x)
    y .= s.region .* x
    return y
end
inputspace(s::SupportProjection) = AnyElementSpace(size(s.region))
outputspace(s::SupportProjection) = AnyElementSpace(size(s.region))
function preallocate_output(
        ::SupportProjection{T, D, A}, x::AbstractArray{S, D}) where {T, S, D, A}
    similar(x)
end

struct FourierTransform{P} <: LinearOperator
    plan::P
    function FourierTransform(dims)
        plan = plan_fft!(zeros(ComplexF64, dims))
        new{typeof(plan)}(plan)
    end
end
struct InvFourierTransform{P} <: LinearOperator
    plan::P
    function InvFourierTransform(dims)
        plan = plan_ifft!(zeros(ComplexF64, dims))
        new{typeof(plan)}(plan)
    end
end
Base.inv(F::FourierTransform) = InvFourierTransform(size(inputspace(F)))
Base.inv(F::InvFourierTransform) = FourierTransform(size(inputspace(F)))

const FT = Union{FourierTransform, InvFourierTransform}
function LinearAlgebra.mul!(y, F::FT, x)
    y .= x
    mul!(y, F.plan, y)
    return y
end
inputspace(F::FT) = AnyElementSpace(size(F.plan))
outputspace(F::FT) = TensorSpace{ComplexF64, length(size(F.plan))}(size(F.plan))
function preallocate_output(F::FT, x::AbstractArray{T, D}) where {T, D}
    @assert size(x) == size(F.plan)
    Array{ComplexF64, D}(undef, size(x))
end
struct Vec{T, D} <: LinearOperator
    m::NTuple{D, Int}
    Vec(m::NTuple{D, Int}, ::Val{T}) where {T, D} = new{T, D}(m)
end
function LinearAlgebra.mul!(
        y::AbstractVector{T},
        v::Vec{T, D},
        x::AbstractArray{T, D}
) where {T, D}
    @assert size(y) == (prod(v.m),)
    fast_reshape!(y, x)
end
inputspace(v::Vec{T, D}) where {T, D} = TensorSpace{T, D}(v.m)
outputspace(v::Vec{T, D}) where {T, D} = TensorSpace{T, 1}((prod(v.m),))
function preallocate_output(v::Vec{T, D}, x::AbstractArray{T, D}) where {T, D}
    similar(x, T, (prod(v.m),))
end
preallocate_input(v::Vec{T, D}) where {T, D} = Array{T, D}(undef, v.m)

struct invVec{T, D} <: LinearOperator
    m::NTuple{D, Int}
end
Base.inv(v::Vec{T, D}) where {T, D} = invVec{T, D}(v.m)
function LinearAlgebra.mul!(
        y::AbstractArray{T, D},
        v::invVec{T, D},
        x::AbstractVector{T}
) where {T, D}
    @assert size(y) == v.m
    fast_reshape!(y, x)
end
function fast_reshape!(y, x)
    @assert length(y) == length(x)
    @inbounds for i in eachindex(x, y)
        y[i] = x[i]
    end
    return y
end

inputspace(v::invVec{T, D}) where {T, D} = TensorSpace{T, 1}((prod(v.m),))
outputspace(v::invVec{T, D}) where {T, D} = TensorSpace{T, D}(v.m)
function preallocate_output(v::invVec{T, D}, x::AbstractArray{T, 1}) where {T, D}
    similar(x, T, v.m)
end
preallocate_input(v::invVec{T, D}) where {T, D} = Array{T, 1}(undef, (prod(v.m),))

struct IdentityOperator{A} <: LinearOperator
    inputspace::A
end
inputspace(id::IdentityOperator) = id.inputspace
outputspace(id::IdentityOperator) = id.inputspace
function LinearAlgebra.mul!(y::AbstractArray, id::IdentityOperator, x::AbstractArray)
    @assert size(x) == size(y) == size(inputspace(id))
    y .= x
    return y
end
preallocate_output(::IdentityOperator, x::AbstractArray) = similar(x)
function preallocate_input(v::IdentityOperator{TensorSpace{T, D}}) where {T, D}
    Array{T, D}(undef, size(inputspace(v)))
end

struct RealPart{T <: Real} <: LinearOperator end
function LinearAlgebra.mul!(y, ::RealPart, x)
    y .= real.(x)
    return y
end

inputspace(::RealPart) = AnySpace()
outputspace(::RealPart{T}) where {T} = AnySizeSpace{T}()
function preallocate_output(
        ::RealPart{T},
        x::AbstractArray{S, D}
) where {T, D, S <: Union{Complex{T}, T}}
    Array{real(T), D}(undef, size(x))
end
