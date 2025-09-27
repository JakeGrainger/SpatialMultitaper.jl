"""
Trait system for handling different process structures in estimates.

- SingleProcessTrait: One process → D-dimensional array of numbers
- MultipleVectorTrait: Multiple processes as vector → (D+2)-dimensional array
- MultipleTupleTrait: Multiple processes as tuple → D-dimensional array of SMatrix
"""
abstract type ProcessTrait end
struct SingleProcessTrait <: ProcessTrait end
struct MultipleVectorTrait <: ProcessTrait end
struct MultipleTupleTrait <: ProcessTrait end

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
    process_trait(data::MultipleSpatialDataVec)

Determine the trait for input spatial data.
"""
process_trait(::MultipleSpatialDataVec) = MultipleVectorTrait()
process_trait(::MultipleSpatialDataTuple) = MultipleTupleTrait()
process_trait(::PointPattern) = SingleProcessTrait()
process_trait(::MarkedPointPattern) = SingleProcessTrait()
process_trait(::GriddedData) = SingleProcessTrait()

function index_process_trait(estimate, i, j)
    index_process_trait(process_trait(estimate), i, j)
end

# single index always is single process
index_process_trait(::MultipleVectorTrait, ::Int, ::Int) = SingleProcessTrait()
index_process_trait(::MultipleTupleTrait, ::Int, ::Int) = SingleProcessTrait()
index_process_trait(::SingleProcessTrait, ::Int, ::Int) = SingleProcessTrait()

# other kinds of indexing throw errors
function index_process_trait(::SingleProcessTrait, i, j)
    throw(ArgumentError("Single process trait only supports single index, got ($i, $j)"))
end

# otherwise returns the same trait as indexing the whole object
index_process_trait(::MultipleVectorTrait, i, j) = MultipleVectorTrait()
index_process_trait(::MultipleTupleTrait, i, j) = MultipleTupleTrait()
