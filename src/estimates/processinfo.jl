struct ProcessInformation{D, T, I1, I2, M, A}
    process_indices_1::I1
    process_indices_2::I2
    mean_product::M
    atoms::A
    function ProcessInformation{D, T}(process_indices_1::I1, process_indices_2::I2,
            mean_product::M, atoms::A) where {D, T, I1, I2, M, A}
        @argcheck length(process_indices_1) == size(mean_product, 1)
        @argcheck length(process_indices_2) == size(mean_product, 2)
        @argcheck size(mean_product) == size(atoms)
        new{D, T, I1, I2, M, A}(process_indices_1, process_indices_2, mean_product, atoms)
    end
    function ProcessInformation(data::SpatialData; mean_method = DefaultMean())
        process_indices_1 = copy(propertynames(data))
        process_indices_2 = copy(propertynames(data))
        λ = mean_estimate(data, mean_method)
        mean_product = λ * λ'
        zero_atom = covariance_zero_atom(data)
        D = embeddim(data)
        trait = process_trait(data)
        ProcessInformation{D, typeof(trait)}(
            process_indices_1, process_indices_2, mean_product, zero_atom)
    end
end
Meshes.embeddim(::ProcessInformation{D}) where {D} = D

function _process_trait_from_info(::ProcessInformation{D, SingleProcessTrait}) where {D}
    SingleProcessTrait()
end
function _process_trait_from_info(::ProcessInformation{D, MultipleVectorTrait}) where {D}
    MultipleVectorTrait()
end
function _process_trait_from_info(::ProcessInformation{D, MultipleTupleTrait}) where {D}
    MultipleTupleTrait()
end

"""
    checkprocessinformation(processinformation, est)

Validate that the ProcessInformation is consistent with the estimate array dimensions.

Returns `(P, Q)` where:
- `P` is the number of processes in the first set
- `Q` is the number of processes in the second set

# Logic
1. Determines process counts from the structure of process_indices_1 and process_indices_2
2. Validates estimate array dimensions against expected structure
3. For SMatrix estimates, validates against matrix dimensions

# Process Count Rules
- Vector/Tuple indices → multiple processes (count = length)
- Single indices → single process (count = 1)
"""
function checkprocessinformation(processinformation::ProcessInformation, est)
    # Determine process counts from index structure
    P = _count_processes(processinformation.process_indices_1)
    Q = _count_processes(processinformation.process_indices_2)

    # Validate dimensions based on estimate type and process configuration
    _validate_estimate_dimensions(
        processinformation, est, P, Q, _process_trait_from_info(processinformation))

    return P, Q
end

"""
    _count_processes(indices)

Count the number of processes from the index structure.
"""
_count_processes(indices::Union{AbstractVector, Tuple}) = length(indices)
_count_processes(indices) = 1  # Single process case

"""
    _validate_estimate_dimensions(processinformation, est, P, Q)

Validate that estimate dimensions are consistent with process configuration.
"""
function _validate_estimate_dimensions(
        processinformation::ProcessInformation{D},
        est::AbstractArray{<:SMatrix{P_mat, Q_mat}},
        P, Q, ::MultipleTupleTrait
) where {D, P_mat, Q_mat}
    # SMatrix case: validate against matrix dimensions
    if P != P_mat
        throw(ArgumentError("Process count mismatch: expected $P processes in first set, got SMatrix with $P_mat rows"))
    end
    if Q != Q_mat
        throw(ArgumentError("Process count mismatch: expected $Q processes in second set, got SMatrix with $Q_mat columns"))
    end

    # Check spatial dimensions
    spatial_dims = ndims(est)
    if spatial_dims != D && spatial_dims != 1
        throw(DimensionMismatch("Expected $D spatial dimensions or 1 (collapsed), got $spatial_dims"))
    end

    return nothing
end

function _validate_estimate_dimensions(
        ::ProcessInformation{D},
        est::AbstractArray{T, N},
        P, Q, ::MultipleVectorTrait
) where {D, T, N}
    if !(ndims(est) == D + 2 || ndims(est) == 3)  # Allow for collapsed spatial dims
        throw(DimensionMismatch("Expected $(D+2) dimensions or 3 (if isotropic), got $(ndims(est))"))
    end
    if size(est, 1) != P
        throw(ArgumentError("First dimension size $(size(est,1)) does not match process count $P"))
    end
    if size(est, 2) != Q
        throw(ArgumentError("Second dimension size $(size(est,2)) does not match process count $Q"))
    end
    return nothing
end

function _validate_estimate_dimensions(
        ::ProcessInformation{D},
        est::AbstractArray{T, N},
        P, Q, ::SingleProcessTrait
) where {D, T, N}
    if !(ndims(est) == D || ndims(est) == 1)  # Allow for collapsed spatial dims
        throw(DimensionMismatch("Expected $(D) dimensions or 1 (if isotropic), got $(ndims(est))"))
    end
    if !(P == Q == 1)
        throw(ArgumentError("Single process should only be indexed at (1,1), but got P=$P, Q=$Q"))
    end
    return nothing
end

function is_same_process_sets(processinformation::ProcessInformation)
    p1 = processinformation.process_indices_1
    p2 = processinformation.process_indices_2
    p1 == p2
end
