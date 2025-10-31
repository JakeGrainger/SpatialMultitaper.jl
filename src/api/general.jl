function rotational_estimate(est::AnisotropicEstimate{E, D}; kwargs...) where {E, D}
    return compute(RotationalEstimate{E, D, typeof(est)}, est; kwargs...)
end
