struct MarginallyTransformedEstimate{E, D, P, Q, N, S, A, T, IP, IE, F} <:
       AbstractEstimate{E, D, P, Q, N}
    argument::A
    estimate::T
    processinformation::IP
    estimationinformation::IE
    function MarginallyTransformedEstimate{E, S, N, F}(argument::A, estimate::T,
            processinfo::ProcessInformation{D}, estimationinfo::IE) where {
            E, S, N, F, A, T, D, IE}
        P, Q = checkinputs(argument, estimate, processinfo)
        IP = typeof(processinfo)
        new{E, D, P, Q, N, S, A, T, IP, IE, F}(
            argument, estimate, processinfo, estimationinfo)
    end
end

function getestimatename(T::Type{<:MarginallyTransformedEstimate})
    "$(gettransformname(T))($(getestimatename(getoriginaltype(T))))"
end

function getestimate(est::MarginallyTransformedEstimate)
    est.estimate
end
function getargument(est::MarginallyTransformedEstimate)
    est.argument
end
function _construct_estimate_subset(
        ::Type{<:MarginallyTransformedEstimate{E, D, P, Q, N, S, A, T, IP, IE, F}},
        trait::Type{<:EstimateTrait},
        args...
) where {E, D, P, Q, N, S, A, T, IP, IE, F}
    return MarginallyTransformedEstimate{trait, S, N, F}(args...)
end

# Helper functions for dispatch on original type
function getoriginaltype(::Type{<:MarginallyTransformedEstimate{
        E, D, P, Q, N, S}}) where {E, D, P, Q, N, S}
    return S
end
gettransformtype(est::MarginallyTransformedEstimate) = gettransformtype(typeof(est))
function gettransformtype(::Type{<:MarginallyTransformedEstimate{
        E, D, P, Q, N, S, A, T, IP, IE, F}}) where {E, D, P, Q, N, S, A, T, IP, IE, F}
    F
end

gettransformname(est::MarginallyTransformedEstimate) = gettransformname(typeof(est))
function gettransformname(T::Type{<:MarginallyTransformedEstimate})
    return lstrip(string(nameof(gettransformtype(T))), '#')
end

function apply_marginal_transform(
        transform::F, est::AbstractEstimate{E, D, P, Q, N}) where {F, E, D, P, Q, N}
    argument = getargument(est)
    estimate = apply_marginal_transform(transform, getestimate(est))
    processinfo = getprocessinformation(est)
    estimationinfo = getestimationinformation(est)
    S = typeof(est)
    return MarginallyTransformedEstimate{E, S, N, F}(
        argument, estimate, processinfo, estimationinfo)
end

function apply_marginal_transform(transform::F, x::AbstractArray{<:SMatrix}) where {F}
    return map(y -> transform.(y), x)
end

function apply_marginal_transform(transform, x::AbstractArray{<:Number})
    transform.(x)
end

Base.real(x::AbstractEstimate) = apply_marginal_transform(real, x)
Base.imag(x::AbstractEstimate) = apply_marginal_transform(imag, x)
Base.conj(x::AbstractEstimate) = apply_marginal_transform(conj, x)
Base.abs(x::AbstractEstimate) = apply_marginal_transform(abs, x)
Base.abs2(x::AbstractEstimate) = apply_marginal_transform(abs2, x)
Base.angle(x::AbstractEstimate) = apply_marginal_transform(angle, x)
Base.log(x::AbstractEstimate) = apply_marginal_transform(log, x)
Base.exp(x::AbstractEstimate) = apply_marginal_transform(exp, x)
