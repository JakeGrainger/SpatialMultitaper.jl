"""
    _dft_to_spectral_matrix(data, J_n, nk)

Compute the spectral matrix from DFTs for different process types.

# Arguments
- `J_n`: The tapered DFTs, with shape depending on the process type:
    - For `SingleProcessTrait`: `n₁ × ... × n_D × M` array
    - For `MultipleVectorTrait`: `P × M × n₁ × ... × n_D` array
    - For `MultipleTupleTrait`: A `n₁ × ... × n_D × M` array of `SVector`s of length `P`
- Trait type: Specifies the process structure.

# Returns
- The spectral matrix or array of spectral matrices, with shape depending on the input.

# Notes
- For single process: returns an array of power spectral densities.
- For multiple processes (tuple): returns an array of spectral matrices.
- For multiple vector processes: returns a D+2 array where each slice is a spectral matrix.
"""
function _dft_to_spectral_matrix!(
        S_mat::AbstractArray, J_n::AbstractArray, ::MultipleVectorTrait)
    for i in CartesianIndices(size(J_n)[3:end])
        @views _compute_spectral_matrix!(S_mat[:, :, i], J_n[:, :, i])
    end
    return S_mat
end

function _dft_to_spectral_matrix!(
        S_mat::AbstractArray, J_n::AbstractArray, ::SingleProcessTrait)
    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(abs2, @view J_n[i, :])
    end
    return S_mat
end
function _dft_to_spectral_matrix!(
        S_mat::Array{<:SMatrix{P, P}}, J_n::AbstractArray{<:SVector{P}, N},
        ::MultipleTupleTrait) where {P, N}
    # Validate dimensions
    @argcheck size(S_mat) == size(J_n)[1:(end - 1)]

    for i in CartesianIndices(S_mat)
        S_mat[i] = mean(_compute_spectral_matrix(J_n[i, m]) for m in axes(J_n, N))
    end
    return S_mat
end

"""
    _compute_spectral_matrix(x::AbstractVector)

Compute the outer product spectral matrix for a vector.
"""
_compute_spectral_matrix(x::AbstractVector) = x * x'

"""
    _compute_spectral_matrix(x::AbstractMatrix)

Compute the averaged spectral matrix for a matrix.
"""
_compute_spectral_matrix(x::SMatrix) = (x * x') ./ size(x, 2)

function _compute_spectral_matrix!(y, x::AbstractMatrix)
    for i in axes(y, 1), j in axes(y, 2)
        y[i, j] = @views mean(x[i, :] .* conj.(x[j, :]))
    end
end

"""
    _make_wavenumber_grid(nk, kmax)

Create a wavenumber grid with dimension-specific parameters.
"""
function _make_wavenumber_grid(nk, kmax)
    wavenumber = _choose_wavenumbers_1d.(nk, kmax)
    return wavenumber
end
