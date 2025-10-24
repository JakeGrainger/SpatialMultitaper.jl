abstract type RotationalKernel end

"""
    GaussKernel{T <: Real}

Gaussian smoothing kernel for rotational averaging.

# Fields
- `bw::T`: Bandwidth parameter controlling the width of the Gaussian

# Mathematical Form
K(x) = exp(-x²/(2σ²)) / σ where σ is the bandwidth.
"""
struct GaussKernel{T <: Real} <: RotationalKernel
    bw::T
end

"""
    (f::GaussKernel)(x)

Evaluate the Gaussian kernel at distance x.
"""
function (f::GaussKernel)(x)
    return exp(-(x^2) / (2 * f.bw^2)) / f.bw
end

"""
    RectKernel{T <: Real}

Rectangular (uniform) smoothing kernel for rotational averaging.

# Fields
- `bw::T`: Bandwidth parameter controlling the width of the rectangular window

# Mathematical Form
K(x) = 1/h if |x/h| < 1/2, 0 otherwise, where h is the bandwidth.
"""
struct RectKernel{T <: Real} <: RotationalKernel
    bw::T
end

"""
    (f::RectKernel)(x)

Evaluate the rectangular kernel at distance x.
"""
function (f::RectKernel)(x)
    return (abs(x / f.bw) < 1 / 2) / f.bw
end
