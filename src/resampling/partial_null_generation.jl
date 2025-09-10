function prediction_kernel(spec::SpectralEstimate)
    kernel_ft = prediction_kernel_ft(spec)
    kernel = prediction_kernel_ft2space.(Ref(kernel_ft.freq), kernel_ft.kernels) # TODO: got here
end

function prediction_kernel_ft2space(freq, power::AbstractArray{D,T}) where {D,T<:Number}
    @assert length(freq) == ndims(power) - 1
    (length(power) * prod(step.(freq))) .*
    fftshift(ifft(ifftshift(power, 2:ndims(power)), 2:ndims(power)), 2:ndims(power))
end

function prediction_kernel_ft2space(freq, power::AbstractArray{S,SMatrix})
    prediction_kernel_ft2space(freq, svectors2array(power))
end

function svectors2array(x::Array{D,SVector{P,T}}) where {D,P,T}
    reinterpret(reshape, T, x)
end


"""
    prediction_kernel_ft(spec)

Computes the Fourier transform of the prediction kernel given estimates of the spectral density function.
"""
function prediction_kernel_ft(spec::SpectralEstimate{D,F,P,N,T}) where {D,F,P,N,T}
    kernels = ntuple(
        idx ->
            apply_transform(single_prediction_kernel_ft, spec.power, idx, spec.ntapers),
        Val{P}(),
    )
    return (freq = spec.freq, kernels = kernels)
end

function single_prediction_kernel_ft(S_mat::AbstractMatrix, idx::Int, ntapers::Int)
    scaling = ntapers / (ntapers - size(S_mat, 1) + 1)
    return single_prediction_kernel_ft(S_mat, idx, nothing) .* scaling
end

"""
    single_prediction_kernel_ft(S_mat::SMatrix{P,P,T,L}, idx::Int, ::Nothing)

Computes the Fourier transform of the prediction kernel for a single variable given an estimate of the spectral density function.
"""
function single_prediction_kernel_ft(
    S_mat::SMatrix{P,P,T,L},
    idx::Int,
    ::Nothing,
) where {P,T,L}
    @boundscheck checkbounds(1:P, idx)
    other_idx =
        StaticArrays.sacollect(SVector{P - 1,Int}, ApplyArray(vcat, 1:idx-1, idx+1:P))
    S_XZ = S_mat[idx, other_idx]
    S_ZZ = S_mat[other_idx, other_idx]
    return S_XZ / S_ZZ
end

"""
    single_prediction_kernel_ft(S_mat::AbstractMatrix, idx::Int, ::Nothing)

Computes the Fourier transform of the prediction kernel for a single variable given an estimate of the spectral density function.
"""
function single_prediction_kernel_ft(S_mat::AbstractMatrix, idx::Int, ::Nothing)
    @boundscheck checkbounds(1:size(S_mat, 1), idx)
    S_XZ = S_mat[idx, Not(idx)]
    S_ZZ = S_mat[Not(idx), Not(idx)]
    return S_XZ / S_ZZ
end