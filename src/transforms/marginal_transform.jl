struct MarginallyTransformedEstimate{A, E, T, S, F, O, D, P, Q, N} <:
       AbstractEstimate{D, P, Q, N}
    argument::A
    estimate::E
    processinformation::T
    estimationinformation::S
    transform::F
    originaltype::O
    function MarginallyTransformedEstimate(
            argument::NTuple{N}, estimate, processinfo::ProcessInformation{D}, estimationinfo,
            transform, originaltype) where {N, D}
        P, Q = checkinputs(argument, estimate, processinfo)
        new{typeof(argument), typeof(estimate),
            typeof(processinfo), typeof(estimationinfo), typeof(transform),
            typeof(originaltype), D, P, Q, N}(
            argument, estimate, processinfo, estimationinfo, transform, originaltype)
    end
end
function getestimate(est::MarginallyTransformedEstimate)
    est.estimate
end
function getargument(est::MarginallyTransformedEstimate)
    est.argument
end
getextrainformation(est::MarginallyTransformedEstimate) = (est.transform, est.originaltype)

function apply_marginal_transform(transform::F, est::AbstractEstimate) where {F}
    MarginallyTransformedEstimate(
        getargument(est), apply_marginal_transform(transform, getestimate(est)),
        getprocessinformation(est), getestimationinformation(est),
        transform, getestimatename(est))
end

function apply_marginal_transform(transform::F, x::AbstractArray{<:SMatrix}) where {F}
    return map(y -> transform.(y), x)
end

function apply_marginal_transform(transform, x::AbstractArray{<:Number})
    transform.(x)
end

Base.real(x::AbstractEstimate) = apply_marginal_transform(real, x)
Base.imag(x::AbstractEstimate) = apply_marginal_transform(imag, x)
Base.abs2(x::AbstractEstimate) = apply_marginal_transform(abs2, x)
Base.abs(x::AbstractEstimate) = apply_marginal_transform(abs, x)
Base.conj(x::AbstractEstimate) = apply_marginal_transform(conj, x)
Base.log(x::AbstractEstimate) = apply_marginal_transform(log, x)
Base.exp(x::AbstractEstimate) = apply_marginal_transform(exp, x)
