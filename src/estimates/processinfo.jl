
struct ProcessInformation{D, T, I1, I2, M, A}
    process_indices_1::I1
    process_indices_2::I2
    mean_product::M
    atoms::A
    function ProcessInformation{D, T}(process_indices_1::I1, process_indices_2::I2,
            mean_product::M, atoms::A) where {D, T, I1, I2, M, A}
        new{D, T, I1, I2, M, A}(process_indices_1, process_indices_2, mean_product, atoms)
    end
    function ProcessInformation(data::SpatialData; mean_method = DefaultMean())
        process_indices_1 = copy(propertynames(data))
        process_indices_2 = copy(propertynames(data))
        λ = mean_estimate(data, mean_method)
        mean_product = λ * λ'
        zero_atom = covariance_zero_atom(data)
        D = embeddim(data)
        trait = estimate_trait(data)
        ProcessInformation{D, typeof(trait)}(
            process_indices_1, process_indices_2, mean_product, zero_atom)
    end
end
Meshes.embeddim(::ProcessInformation{D}) where {D} = D

function checkprocessinformation(
        processinformation::ProcessInformation, P, Q)
    if length(processinformation.process_indices_1) != P
        throw(ArgumentError("processinformation.process_indices_1 should have length $P"))
    end
    if length(processinformation.process_indices_2) != Q
        throw(ArgumentError("processinformation.process_indices_2 should have length $Q"))
    end
    if size(processinformation.mean_product, 1) != P ||
       size(processinformation.mean_product, 2) != Q
        throw(ArgumentError(
            "processinformation.mean_product should have size ($P, $Q), but has size $(size(processinformation.mean_product))."
        ))
    end
    if size(processinformation.atoms, 1) != P || size(processinformation.atoms, 2) != Q
        throw(ArgumentError(
            "processinformation.atoms should have size ($P, $Q), but has size $(size(processinformation.atoms))."
        ))
    end
end

function is_same_process_sets(processinformation::ProcessInformation)
    p1 = processinformation.process_indices_1
    p2 = processinformation.process_indices_2
    p1 == p2
end
