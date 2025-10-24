
"""
    rotational_estimate(est::AnisotropicEstimate{E}; radii, kernel) where {E}

Compute a rotationally averaged estimate from an anisotropic estimate.

Rotational averaging integrates the estimate over circular annuli at different
radii, producing an isotropic result that depends only on distance from the origin.
This is useful for analyzing scale-dependent properties without directional bias.

# Arguments
- `est::AnisotropicEstimate`: The anisotropic estimate to average
- `radii`: Radial distances for evaluation (default: `default_rotational_radii(est)`)
- `kernel`: Smoothing kernel for averaging (default: `default_rotational_kernel(est)`)

# Returns
A `RotationalEstimate` containing the rotationally averaged values.
The exception is if kernel is `NoRotational`, in which case the input estimate is returned.
This is used to skip rotational averaging in some downstream cases.

# Examples
```julia
# Basic rotational averaging
rot_est = rotational_estimate(spectrum)

# With custom parameters
radii = range(0.1, 2.0, length=50)
kernel = GaussKernel(0.1)  # Gaussian smoothing with bandwidth 0.1
rot_est = rotational_estimate(spectrum, radii=radii, kernel=kernel)
```
"""
function rotational_estimate(
        est::AnisotropicEstimate; radii = default_rotational_radii(est),
        kernel = default_rotational_kernel(est))
    return _rotational_estimate(est, radii, kernel)
end
