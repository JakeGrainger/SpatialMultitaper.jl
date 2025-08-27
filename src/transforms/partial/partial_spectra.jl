struct PartialSpectra{D,F,P,N} <: AnisotropicEstimate{D,P}
    freq::NTuple{D,F}
    partial_spectra::N
    function PartialSpectra(freq::NTuple{D,F}, partial_coherence) where {D,F}
        P = checkinputs(freq, partial_coherence)
        new{D,F,P,typeof(partial_coherence)}(freq, partial_coherence)
    end
end
getargument(est::PartialSpectra) = est.freq
getestimate(est::PartialSpectra) = est.partial_spectra


"""
	partial_spectra(x::Matrix)

If only a matrix `x` is passed, computes the partial spectra for all indices. 
The diagonal elements are the spectra of residuals given all other processes
The i j th element (for i ≠ j) is the partial spectra between the ith and jth process given all other processes not including i and j.
	
	partial_spectra(x::Matrix, i1::Int, i2::Int, c1, c2)

If specific indices are requested, computes the partial spectra for the i1th index conditioned on the indices in c1 vs the i2th index conditioned on the indices in c2.
"""
function partial_spectra(x::SMatrix, ::Nothing)
    g = inv(x)
    A = diagm((diag(g)))
    g2 = abs2.(g)
    denom = A * ones(typeof(x)) * A - g2 + diagm(diag(g2))
    return (g ./ denom) .* (2I - ones(typeof(x)))
    # computes -gⱼₖ / (gⱼⱼ gₖₖ - |gⱼₖ|²) if j ≠ k
    # computes  1 / gⱼⱼ if j = k
end

function partial_spectra(x::AbstractMatrix, ::Nothing)
    C = inv(x)

    return [
        i == j ? 1 / C[i, i] : -C[i, j] / (C[i, i] * C[j, j] - abs2(C[i, j])) for
        i in axes(C, 1), j in axes(C, 2)
    ]
end

function partial_spectra(x::SMatrix{2,2,T,4}, ::Nothing) where {T}
    g = inv(x)
    p11 = 1 / g[1, 1]
    p22 = 1 / g[2, 2]
    p12 = -g[1, 2] / (g[1, 1] * g[2, 2] - abs2(g[1, 2]))
    p21 = conj(p12)

    return SMatrix{2,2,T,4}(
        p11,
        p21,
        p12,
        p22, # column major
    )
end

function partial_spectra(x::SMatrix{Q,Q,T,N}, ntapers::Int) where {Q,T,N}
    p = partial_spectra(x, nothing)
    denom = ntapers .* ones(typeof(x)) .- Q .+ 2 - I # so that M - Q + 2 off diag and M - Q + 1 on diag
    return ntapers ./ denom .* p
end

function partial_spectra(x::AbstractMatrix{T}, ntapers::Int) where {T}
    Q = size(x, 1)
    p = partial_spectra(x, nothing)
    for i in axes(p, 1), j in axes(p, 2)
        @inbounds p[i, j] *= ntapers / (ntapers - Q + 2 - (i == j))
    end
    return p
end

partial_spectra(spectrum::SpectralEstimate; partial_type::PartialType = UsualPartial()) =
    _partial_spectra(spectrum, partial_type)
function _partial_spectra(spectrum::SpectralEstimate, ::UsualPartial)
    return PartialSpectra(
        spectrum.freq,
        apply_transform(partial_spectra, spectrum.power, spectrum.ntapers),
    )
end
function _partial_spectra(spectrum::SpectralEstimate, ::SplitPartial)
    return PartialSpectra(
        spectrum.freq,
        apply_transform(split_partial_spectra, spectrum.power, spectrum.ntapers),
    )
end

function partial_spectra_uncorrected(spectrum::SpectralEstimate)
    return PartialSpectra(
        spectrum.freq,
        apply_transform(partial_spectra, spectrum.power, nothing),
    )
end

function partial_spectra_uncorrected(spectrum::SpectralEstimate, i1::Int, i2::Int, c1, c2)
    function wrapped_partial_spectra(x)
        partial_spectra(x, i1, i2, c1, c2, nothing)
    end # maybe better to make a version of apply_transform for passing args...
    return PartialSpectra(
        spectrum.freq,
        apply_transform(wrapped_partial_spectra, spectrum.power),
    )
end

##
# debiasing in general isnt known
# function partial_spectra(x::AbstractMatrix, i1::Int, i2::Int, c1, c2, ::Nothing)
#     xi1c1_rdiv_c1c1 = (transpose(x[i1, c1]) / x[c1, c1])
#     xc2c2_ldiv_c2i2 = (x[c2, c2] \ x[c2, i2])

#     return x[i1, i2] - xi1c1_rdiv_c1c1 * x[c1, i2] -
#            transpose(x[i1, c2]) * xc2c2_ldiv_c2i2 +
#            xi1c1_rdiv_c1c1 * x[c1, c2] * xc2c2_ldiv_c2i2
# end

# split_partial_spectra(x::AbstractMatrix, i1::Int, i2::Int, c1, c2, ::Nothing) =
#     partial_spectra(x, i1, i2, c1, c2, nothing)

# function split_partial_spectra(x::AbstractMatrix, i1::Int, i2::Int, c1, c2, ntapers)
#     return split_partial_spectra(x, i1, i2, c1, c2, nothing) * ntapers /
#            (ntapers - length(c1) - length(c2))
# end

function split_partial_spectra(x::Matrix{T}, ntapers) where {T}
    # assumes that the spectrum is for an even number of processes, and we want to give back all partial spectra of the form 
    # i on 1:P not i j with j + P on P+1:2P not i + P j + P
    P = size(x, 1) ÷ 2
    out = zeros(T, P, P)
    firsthalf = 1:P
    for i = 1:P
        temp = partial_spectra(
            view(x, [firsthalf[Not(i)]; i + P], [firsthalf[Not(i)]; i + P]),
            ntapers,
        )
        out[Not(i), i] = @view temp[end, 1:end-1]
        temp_diag =
            partial_spectra(view(x, [firsthalf; i + P], [firsthalf; i + P]), ntapers)
        out[i, i] = temp_diag[i, end] # shifted version of process always put last
    end
    return out
end

function _split_partial_spectra_static(
    x::SMatrix{Q,Q,T,N},
    ntapers,
    ::Val{P},
) where {Q,T,N,P}
    reduce(
        hcat,
        ntuple(
            i -> _split_partial_spectra_static_single_row(x, ntapers, Val{P}(), i),
            Val{P}(),
        ),
    )
end

function _split_partial_spectra_static_single_row(
    x::SMatrix{Q,Q,T,N},
    ntapers,
    ::Val{P},
    i,
) where {Q,T,N,P}
    temp_idx = StaticArrays.sacollect(SVector{P,Int}, ApplyArray(vcat, 1:i-1, i+1:P, P + i))
    temp = partial_spectra(x[temp_idx, temp_idx], ntapers)
    temp_diag_idx = StaticArrays.sacollect(SVector{P + 1,Int}, ApplyArray(vcat, 1:P, P + i))
    temp_diag = partial_spectra(x[temp_diag_idx, temp_diag_idx], ntapers)
    out = [temp[end, 1:end-1]; temp_diag[i, end]]
    idx = StaticArrays.sacollect(SVector{P,Int}, ApplyArray(vcat, 1:i-1, P, i:P-1)) # because we want to put 1:i-1 from temp, then i from temp_diag, then the remainder of temp (i to P-1)
    return out[idx]
end

##

# function split_partial_spectra(x::SMatrix{4,4,T,16}, ntapers) where {T}
#     # assumes that the spectrum is for an even number of processes, and we want to give back all partial spectra of the form 
#     # i on 1:P not i j with j + P on P+1:2P not i + P j + P


#     p11 = split_partial_spectra(x, 1, 3, 2, 4, ntapers)
#     p12 = x[1, 4]
#     p22 = split_partial_spectra(x, 2, 4, 1, 3, ntapers)

#     p21 = conj(p12)
#     return SMatrix{2,2,T,4}(p11, p21, p12, p22)
# end

# function split_partial_spectra(x::SMatrix{6,6,T,36}, ntapers) where {T}
#     # assumes that the spectrum is for an even number of processes, and we want to give back all partial spectra of the form 
#     # i on 1:P not i j with j + P on P+1:2P not i + P j + P

#     p11 = split_partial_spectra(x, 1, 4, SVector(2, 3), SVector(5, 6), ntapers)
#     p12 = split_partial_spectra(x, 1, 5, 3, 6, ntapers)
#     p13 = split_partial_spectra(x, 1, 6, 2, 5, ntapers)
#     p22 = split_partial_spectra(x, 2, 5, SVector(1, 3), SVector(4, 6), ntapers)
#     p23 = split_partial_spectra(x, 2, 6, 1, 4, ntapers)
#     p33 = split_partial_spectra(x, 3, 6, SVector(1, 2), SVector(4, 5), ntapers)

#     p21 = conj(p12)
#     p31 = conj(p13)
#     p32 = conj(p23)
#     return SMatrix{3,3,T,9}(p11, p21, p31, p12, p22, p32, p13, p23, p33)
# end

# function split_partial_spectra(x::SMatrix{8,8,T,64}, ntapers) where {T}
#     # assumes that the spectrum is for an even number of processes, and we want to give back all partial spectra of the form 
#     # i on 1:P not i j with j + P on P+1:2P not i + P j + P

#     p11 = split_partial_spectra(x, 1, 5, SVector(2, 3, 4), SVector(6, 7, 8), ntapers)
#     p12 = split_partial_spectra(x, 1, 6, SVector(3, 4), SVector(7, 8), ntapers)
#     p13 = split_partial_spectra(x, 1, 7, SVector(2, 4), SVector(6, 8), ntapers)
#     p14 = split_partial_spectra(x, 1, 8, SVector(2, 3), SVector(6, 7), ntapers)
#     p22 = split_partial_spectra(x, 2, 6, SVector(1, 3, 4), SVector(5, 7, 8), ntapers)
#     p23 = split_partial_spectra(x, 2, 7, SVector(1, 4), SVector(5, 8), ntapers)
#     p24 = split_partial_spectra(x, 2, 8, SVector(1, 3), SVector(5, 7), ntapers)
#     p33 = split_partial_spectra(x, 3, 7, SVector(1, 2, 4), SVector(5, 6, 8), ntapers)
#     p34 = split_partial_spectra(x, 3, 8, SVector(1, 2), SVector(5, 6), ntapers)
#     p44 = split_partial_spectra(x, 4, 8, SVector(1, 2, 3), SVector(5, 6, 7), ntapers)

#     p21 = conj(p12)
#     p31 = conj(p13)
#     p41 = conj(p14)
#     p32 = conj(p23)
#     p42 = conj(p24)
#     p43 = conj(p34)

#     return SMatrix{4,4,T}(
#         p11,
#         p21,
#         p31,
#         p41,
#         p12,
#         p22,
#         p32,
#         p42,
#         p13,
#         p23,
#         p33,
#         p43,
#         p14,
#         p24,
#         p34,
#         p44,
#     )
# end

# split_partial_spectra(x::SMatrix{10,10,T,100}, ntapers) where {T} =
#     _split_partial_spectra_static(x, ntapers, Val{5}())

function split_partial_spectra(x::SMatrix{Q,Q,T,N}, ntapers) where {Q,T,N}
    P = Q ÷ 2
    _split_partial_spectra_static(x, ntapers, Val{P}())
end

# function _split_partial_spectra_static(
#     x::SMatrix{Q,Q,T,N},
#     ntapers,
#     ::Val{P},
# ) where {Q,T,N,P}
#     function x_to_partial(i, j)
#         split_partial_spectra(
#             x,
#             i,
#             P + j,
#             view(1:P, Not(SVector(i, j))),
#             view((P+1):2P, Not(SVector(i, j))),
#             ntapers,
#         )
#     end
#     x_to_partial.(transpose(SOneTo(P)), SOneTo(P))
# end

function split_partial_spectra(spectrum::SpectralEstimate)
    return PartialSpectra(
        spectrum.freq,
        apply_transform(split_partial_spectra, spectrum.power, spectrum.ntapers),
    )
end