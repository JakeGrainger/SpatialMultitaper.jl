"""
    is_same_process_sets(est::AbstractEstimate) -> Bool

Check if an estimate uses the same process sets by delegating to the
process information.

# Arguments
- `est::AbstractEstimate`: The estimate object

# Returns
- `Bool`: Whether the estimate uses same process sets
"""
function is_same_process_sets(est::AbstractEstimate)
    is_same_process_sets(get_process_information(est))
end

function Base.:(==)(a::AbstractEstimate, b::AbstractEstimate)
    get_estimates(a) == get_estimates(b) &&
        get_evaluation_points(a) == get_evaluation_points(b) &&
        get_process_information(a) == get_process_information(b) &&
        get_estimation_information(a) == get_estimation_information(b) &&
        get_extra_information(a) == get_extra_information(b)
end

Meshes.embeddim(::AbstractEstimate{E, D}) where {E, D} = D
Base.size(est::AbstractEstimate) = size(get_process_information(est).mean_product)
function processnames(estimate::AbstractEstimate)
    processinfo = get_process_information(estimate)
    (processinfo.process_indices_1, processinfo.process_indices_2)
end

"""
    process_trait(estimate::AbstractEstimate)

Determine the trait for an estimate based on its process information.
"""
function process_trait(estimate::AbstractEstimate)
    process_info = get_process_information(estimate)
    _process_trait_from_info(process_info)
end

function getprocessinformationindex(est::AbstractEstimate, p, q)
    processinformation = get_process_information(est)
    processtrait = index_process_trait(est, p, q)
    return _getprocessinformationindex(processinformation, processtrait, p, q)
end

function _getprocessinformationindex(
        processinformation, processtrait::SingleProcessTrait, p, q)
    return ProcessInformation{embeddim(processinformation), typeof(processtrait)}(
        processinformation.process_indices_1[p],
        processinformation.process_indices_2[q],
        processinformation.mean_product[p, q],
        processinformation.atoms[p, q]
    )
end
function _getprocessinformationindex(
        processinformation, ::MultipleVectorTrait, p::Int, q::Int)
    _getprocessinformationindex(processinformation, SingleProcessTrait(), p, q)
end
function _getprocessinformationindex(
        processinformation, ::MultipleTupleTrait, p::Int, q::Int)
    _getprocessinformationindex(processinformation, SingleProcessTrait(), p, q)
end
function _getprocessinformationindex(
        processinformation, processtrait::MultipleVectorTrait, p, q)
    return ProcessInformation{embeddim(processinformation), typeof(processtrait)}(
        processinformation.process_indices_1[p],
        processinformation.process_indices_2[q],
        processinformation.mean_product[_index_to_vec(p), _index_to_vec(q)],
        processinformation.atoms[_index_to_vec(p), _index_to_vec(q)]
    )
end

function _getprocessinformationindex(
        processinformation, processtrait::MultipleTupleTrait, p, q)
    return ProcessInformation{embeddim(processinformation), typeof(processtrait)}(
        processinformation.process_indices_1[p],
        processinformation.process_indices_2[q],
        processinformation.mean_product[_index_to_svec(p), _index_to_svec(q)],
        processinformation.atoms[_index_to_svec(p), _index_to_svec(q)]
    )
end

struct EstimationInformation{T <: Union{Nothing, Int}}
    ntapers::T
end

"""
    collect(estimate::AbstractEstimate)

Produces a Tuple containing the arguments and the estimate.
"""
function Base.collect(estimate::AbstractEstimate)
    tuple(get_evaluation_points_tuple(estimate)..., get_estimates(estimate))
end

function get_evaluation_points_tuple(estimate::AbstractEstimate)
    get_evaluation_points(estimate)::Tuple
end

function get_evaluation_points_tuple(estimate::IsotropicEstimate)
    _argument_to_tuple(get_evaluation_points(estimate))
end
_argument_to_tuple(x::Tuple) = x
_argument_to_tuple(x) = (x,)

## constructor input checking
function checkinputs(
        argument::AbstractVector, estimate, processinformation::ProcessInformation)
    checkinputs((argument,), estimate, processinformation)
end

function checkinputs(argument::NTuple{D}, estimate::AbstractArray{T, N},
        processinformation::ProcessInformation) where {D, T, N}
    @argcheck length(argument) <= ndims(estimate) <= length(argument) + 2
    @argcheck length.(argument) == size(estimate)[(N + 1 - D):end]
    return checkprocessinformation(processinformation, estimate)
end

function checkinputs(argument::Union{NTuple{N, <:Number}, Number}, estimate::Number,
        processinformation::ProcessInformation) where {N}
    return nothing
end
