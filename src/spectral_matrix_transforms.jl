function _transform_spectral_estimate(transform::T, mt_est::SpectralEstimate{F, P, Nothing}) where {T, F, P}
	transformed_power = mapslices(transform, mt_est.power, dims = (1, 2))
	return (freq=mt_est.freq, transformed_power=transformed_power, transformed_power_jackknifed=nothing)
end

function _transform_spectral_estimate(transform::T, mt_est::SpectralEstimate{F, P, Vector{P}}) where {T, F, P}
	transformed_power = mapslices(transform, mt_est.power, dims = (1, 2))
	transformed_power_jackknifed = mapslices.(Ref(transform), mt_est.power_jackknifed, dims = (1, 2))
	return (freq=mt_est.freq, transformed_power=transformed_power, transformed_power_jackknifed=transformed_power_jackknifed)
end
function transform_spectral_estimate(transform, mt_est::SpectralEstimate, args...; kwargs...)
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
		esc(quote
			$x
			function $(x.args[1].args[1])($(x.args[1].args[2]), mt_est::SpectralEstimate, $(x.args[1].args[4:end]...))
				transformed_estimate = transform_spectral_estimate($(process_spectraltransform_kwargs(x.args[1].args[2])), $(x.args[1].args[1]), mt_est, $(process_spectraltransform_arg.(x.args[1].args[4:end])...))
				return (freq=transformed_estimate.freq, $(x.args[1].args[1])=transformed_estimate.transformed_power, $(Symbol("$(x.args[1].args[1])_jackknifed"))=transformed_estimate.transformed_power_jackknifed)
			end
		end)
	else
		esc(quote
			$x
			function $(x.args[1].args[1])(mt_est::SpectralEstimate, $(x.args[1].args[3:end]...))
				transformed_estimate = transform_spectral_estimate($(x.args[1].args[1]), mt_est, $(process_spectraltransform_arg.(x.args[1].args[3:end])...))
				return (freq=transformed_estimate.freq, $(x.args[1].args[1])=transformed_estimate.transformed_power, $(Symbol("$(x.args[1].args[1])_jackknifed"))=transformed_estimate.transformed_power_jackknifed)
			end
		end)
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

@spectraltransform function complex_coherence(x::Matrix)
	return [x[i, j] / sqrt(x[i, i] * x[j, j]) for i ∈ axes(x, 1), j ∈ axes(x, 2)]
end

@spectraltransform function magnitude_coherence(x::Matrix)
	return abs.(complex_coherence(x))
end

@spectraltransform function magnitude_sq_coherence(x::Matrix)
	return abs2.(complex_coherence(x))
end

@spectraltransform function group_delay(x::Matrix)
	return angle.(complex_coherence(x))
end

@spectraltransform function partial_complex_coherence(x::Matrix)
	return -complex_coherence(inv(x))
end

@spectraltransform function partial_magnitude_coherence(x::Matrix)
	return abs.(partial_complex_coherence(x))
end

@spectraltransform function partial_magnitude_sq_coherence(x::Matrix)
	return abs2.(partial_complex_coherence(x))
end

@spectraltransform function partial_group_delay(x::Matrix)
	return angle.(partial_complex_coherence(x))
end

@spectraltransform function partial_spectra(x::Matrix)
	C = inv(x)
	par_coh = -complex_coherence(C)
	Sa = inv.(diag(C))
	return [
		i == j ? Sa[i] : par_coh[i, j] / (1 - abs2(par_coh[i, j])) * sqrt(Sa[i] * Sa[j])
		for i in axes(C, 1), j in axes(C, 2)]
end

@spectraltransform function partial_spectra(x::Matrix, i1::Int, i2::Int, c1, c2) # note, this can be made much more memory efficient
	x[i1,i2] - transpose(x[i1, c1]) / x[c1, c1] * x[c1, i2] - transpose(x[i1, c2]) / x[c2, c2] * x[c2, i2] + transpose(x[i1, c1]) / x[c1, c1] * x[c1, c2] / x[c2, c2] * x[c2, i2]
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