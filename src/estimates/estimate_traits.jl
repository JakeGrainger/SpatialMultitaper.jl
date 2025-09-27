"""
Trait system for handling different process structures in estimates.

- SingleProcessTrait: One process → D-dimensional array of numbers
- MultipleVectorTrait: Multiple processes as vector → (D+2)-dimensional array
- MultipleTupleTrait: Multiple processes as tuple → D-dimensional array of SMatrix
"""
abstract type EstimateTrait end
struct SingleProcessTrait <: EstimateTrait end
struct MultipleVectorTrait <: EstimateTrait end
struct MultipleTupleTrait <: EstimateTrait end

function _estimate_trait_from_info(::ProcessInformation{D, SingleProcessTrait}) where {D}
    SingleProcessTrait()
end
function _estimate_trait_from_info(::ProcessInformation{D, MultipleVectorTrait}) where {D}
    MultipleVectorTrait()
end
function _estimate_trait_from_info(::ProcessInformation{D, MultipleTupleTrait}) where {D}
    MultipleTupleTrait()
end

"""
    estimate_trait(data::MultipleSpatialDataVec)

Determine the trait for input spatial data.
"""
estimate_trait(::MultipleSpatialDataVec) = MultipleVectorTrait()
estimate_trait(::MultipleSpatialDataTuple) = MultipleTupleTrait()
estimate_trait(::PointPattern) = SingleProcessTrait()
estimate_trait(::MarkedPointPattern) = SingleProcessTrait()
estimate_trait(::GriddedData) = SingleProcessTrait()

function index_estimation_trait(estimate, i, j)
    index_estimation_trait(estimate_trait(estimate), i, j)
end

# single index always is single process
index_estimation_trait(::MultipleVectorTrait, ::Int, ::Int) = SingleProcessTrait()
index_estimation_trait(::MultipleTupleTrait, ::Int, ::Int) = SingleProcessTrait()
index_estimation_trait(::SingleProcessTrait, ::Int, ::Int) = SingleProcessTrait()

# other kinds of indexing throw errors
function index_estimation_trait(::SingleProcessTrait, i, j)
    throw(ArgumentError("Single process trait only supports single index, got ($i, $j)"))
end

# otherwise returns the same trait as indexing the whole object
index_estimation_trait(::MultipleVectorTrait, i, j) = MultipleVectorTrait()
index_estimation_trait(::MultipleTupleTrait, i, j) = MultipleTupleTrait()
