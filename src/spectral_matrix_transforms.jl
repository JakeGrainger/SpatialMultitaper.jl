function apply_transform(transform, power, ::Val{P}) where {P}
    # static_power = ssquarematrixscopy(power, Val{P}())
    # return postprocess_transform(transform.(static_power))
    return transform.(power)
end
function apply_transform(transform, power, ::Val{nothing})
    mapslices(transform, power, dims = (1, 2))
end
postprocess_transform(x::Array{SMatrix{N,N,T,L},D}) where {T,D,N,L} =
    inv_ssquarematrixscopy(x)
postprocess_transform(x::Array) = x

function ssquarematrixscopy(x::Array{T,D}, ::Val{N}) where {T,D,N}
    size(x, 1) == N || error("sizes mismatch")
    size(x, 2) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    copy(reinterpret(reshape, SMatrix{N,N,T,N^2}, reshape(x, (N * N, size(x)[3:end]...))))
end

function inv_ssquarematrixscopy(x::Array{SMatrix{N,N,T,L},D}) where {T,D,N,L}
    isbitstype(T) || error("use for bitstypes only")
    reshape(copy(reinterpret(T, x)), (N, N, size(x)...))
end

function _transform_spectral_estimate(
    transform::T,
    mt_est::SpectralEstimate{P,F,N,Nothing},
) where {T,P,F,N}
    transformed_power = apply_transform(transform, mt_est.power, Val{P}())
    return (
        freq = mt_est.freq,
        transformed_power = transformed_power,
        transformed_power_jackknifed = nothing,
    )
end

function _transform_spectral_estimate(
    transform::T,
    mt_est::SpectralEstimate{P,F,N,Vector{N}},
) where {T,P,F,N}
    transformed_power = apply_transform(transform, mt_est.power, Val{P}())
    transformed_power_jackknifed =
        apply_transform.(Ref(transform), mt_est.power_jackknifed, Val{P}())
    return (
        freq = mt_est.freq,
        transformed_power = transformed_power,
        transformed_power_jackknifed = transformed_power_jackknifed,
    )
end
function transform_spectral_estimate(
    transform,
    mt_est::SpectralEstimate,
    args...;
    kwargs...,
)
    wrapped_transform(x) = transform(x, args...; kwargs...)
    _transform_spectral_estimate(wrapped_transform, mt_est)
end

"""
	@spectraltransform function f(x::Matrix, args...; kwargs...)
		...
	end

A macro which creates a user defined function `f`, but also a function which applied `f` to a `SpectralEstimate` with the signature
	`f(mt_est::SpectralEstimate, args...; kwargs...)`
"""
macro spectraltransform(x)
    if x.args[1].args[2] isa Expr && x.args[1].args[2].head == :parameters
        esc(
            quote
                $x
                function $(x.args[1].args[1])(
                    $(x.args[1].args[2]),
                    mt_est::SpectralEstimate,
                    $(x.args[1].args[4:end]...),
                )
                    transformed_estimate = transform_spectral_estimate(
                        $(process_spectraltransform_kwargs(x.args[1].args[2])),
                        $(x.args[1].args[1]),
                        mt_est,
                        $(process_spectraltransform_arg.(x.args[1].args[4:end])...),
                    )
                    return (
                        freq = transformed_estimate.freq,
                        $(x.args[1].args[1]) = transformed_estimate.transformed_power,
                        $(Symbol("$(x.args[1].args[1])_jackknifed")) = transformed_estimate.transformed_power_jackknifed,
                    )
                end
            end,
        )
    else
        esc(
            quote
                $x
                function $(x.args[1].args[1])(
                    mt_est::SpectralEstimate,
                    $(x.args[1].args[3:end]...),
                )
                    transformed_estimate = transform_spectral_estimate(
                        $(x.args[1].args[1]),
                        mt_est,
                        $(process_spectraltransform_arg.(x.args[1].args[3:end])...),
                    )
                    return (
                        freq = transformed_estimate.freq,
                        $(x.args[1].args[1]) = transformed_estimate.transformed_power,
                        $(Symbol("$(x.args[1].args[1])_jackknifed")) = transformed_estimate.transformed_power_jackknifed,
                    )
                end
            end,
        )
    end
end
function process_spectraltransform_kwargs(ex::Expr)
    @assert ex.head == :parameters
    ex_out = deepcopy(ex)
    for i in eachindex(ex_out.args)
        if ex_out.args[i] isa Symbol
            ex_out.args[i] = Expr(:kw, ex_out.args[i], ex_out.args[i])
        elseif ex_out.args[i] isa Expr
            if ex_out.args[i].head == :kw
                if ex_out.args[i].args[1] isa Expr && ex_out.args[i].args[1].head == :(::)
                    ex_out.args[i].args[1] = ex_out.args[i].args[1].args[1]
                end
                ex_out.args[i].args[2] = ex_out.args[i].args[1]
            elseif ex_out.args[i].head == :(::)
                ex_out.args[i] = ex_out.args[i].args[1]
            else
                error("Failed to process keyword argument when creating functions")
            end
        else
            error("Failed to process keyword argument when creating functions")
        end
    end
    return ex_out
end
process_spectraltransform_arg(arg::Symbol) = arg
function process_spectraltransform_arg(arg::Expr)
    if arg.head == :kw
        return process_test_arg(arg.args[1])
    elseif arg.head == :(::)
        return arg.args[1]
    else
        return arg
    end
end
##

@spectraltransform function complex_coherence(x::AbstractMatrix)
    return [x[i, j] / sqrt(x[i, i] * x[j, j]) for i ∈ axes(x, 1), j ∈ axes(x, 2)]
end

@spectraltransform function magnitude_coherence(x::AbstractMatrix)
    return abs.(complex_coherence(x))
end

@spectraltransform function magnitude_sq_coherence(x::AbstractMatrix)
    return abs2.(complex_coherence(x))
end

@spectraltransform function phase(x::AbstractMatrix)
    return angle.(complex_coherence(x))
end

@spectraltransform function partial_complex_coherence(x::AbstractMatrix)
    return -complex_coherence(inv(x))
end

@spectraltransform function partial_magnitude_coherence(x::AbstractMatrix)
    return abs.(partial_complex_coherence(x))
end

@spectraltransform function partial_magnitude_sq_coherence(x::AbstractMatrix)
    return abs2.(partial_complex_coherence(x))
end

@spectraltransform function partial_phase(x::AbstractMatrix)
    return angle.(partial_complex_coherence(x))
end

@spectraltransform function partial_spectra(x::AbstractMatrix)
    C = inv(x)
    par_coh = -complex_coherence(C)
    invCd = inv.(diag(C))
    invCd_root = diagm(sqrt.(invCd))
    part1 = (invCd_root * par_coh * invCd_root)
    (part1 - diagm(diag(part1)) + diagm(invCd)) ./ (1.0 .- abs2.(par_coh) + I)
end


@spectraltransform function partial_spectra(x::AbstractMatrix, i1::Int, i2::Int, c1, c2) # note, this can be made much more memory efficient
    x[i1, i2] - transpose(x[i1, c1]) / x[c1, c1] * x[c1, i2] -
    transpose(x[i1, c2]) / x[c2, c2] * x[c2, i2] +
    transpose(x[i1, c1]) / x[c1, c1] * x[c1, c2] / x[c2, c2] * x[c2, i2]
end

"""
	partial_spectra(x::Matrix)

If only a matrix `x` is passed, computes the partial spectra for all indices. 
The diagonal elements are the spectra of residuals given all other processes
The i j th element (for i ≠ j) is the partial coherence between the ith and jth process given all other processes not including i and j.
	
	partial_spectra(x::Matrix, i1::Int, i2::Int, c1, c2)

If specific indices are requested, computes the partial spectra for the i1th index conditioned on the indices in c1 vs the i2th index conditioned on the indices in c2.
"""
partial_spectra
