struct PartialKFunction{R,T,D,P} <: IsotropicEstimate{D,P}
    radii::R
    partial_K_function::T
    function PartialKFunction(radii::R, K::T, ::Val{D}) where {R,T,D}
        P = checkinputs(radii, K)
        new{R,T,D,P}(radii, K)
    end
end

getargument(f::PartialKFunction) = f.radii
getestimate(f::PartialKFunction) = f.partial_K_function
getextrafields(::PartialKFunction{R,T,D,P}) where {R,T,D,P} = (Val{D}(),)

function partial_K_function(
    c::PartialCFunction{R,T,D,1},
    λ;
    partial_type = UsualPartial(),
) where {R,T,D}
    invλ = 1 / λ[1]
    return PartialKFunction(
        c.radii,
        C2K(c.radii, c.partial_C_function, invλ, invλ, Val{D}()),
        Val{D}(),
    )
end

function partial_K_function(
    c::PartialCFunction{R,T,D,P},
    λ::NTuple{Q,<:Number};
    partial_type = UsualPartial(),
) where {R,T<:AbstractArray,D,P,Q}
    if partial_type isa UsualPartial
        @assert Q == P "C function and λ should be the same size for UsualPartial"
    else
        @assert Q == P * 2 "λ should be twice the size of C function for  SplitPartial"
    end

    invλ1, invλ2 = prepare_inverse_lambda_for_K_function(λ, partial_type)
    return PartialKFunction(
        c.radii,
        C2K(c.radii, c.partial_C_function, invλ1, invλ2, Val{D}()),
        Val{D}(),
    )
end

function prepare_inverse_lambda_for_K_function(λ, ::UsualPartial)
    invλ = diagm(1 ./ SVector(λ...))
    return invλ, invλ
end

function prepare_inverse_lambda_for_K_function(λ, ::SplitPartial)
    invλ1 = diagm(1 ./ SVector(λ[1:end÷2]...))
    invλ2 = diagm(1 ./ SVector(λ[end÷2+1:end]...))
    return invλ1, invλ2
end

function partial_K_function(
    data,
    region;
    radii,
    nfreq,
    fmax,
    tapers,
    mean_method::MeanEstimationMethod = DefaultMean(),
    partial_type::PartialType = UsualPartial(),
)
    c = partial_C_function(
        data,
        region;
        radii = radii,
        nfreq = nfreq,
        fmax = fmax,
        tapers = tapers,
        mean_method = mean_method,
        partial_type = partial_type,
    )
    λ = mean_estimate(data, region, mean_method)
    return partial_K_function(c, λ, partial_type = partial_type)
end