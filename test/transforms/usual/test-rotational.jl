mt_est = make_simple_example()
radii = 0:0.1:2
kernel = Spmt.RectKernel(0.3)
rot_est = rotational_estimate(mt_est, radii = radii, kernel = kernel)
@test rot_est isa Spmt.RotationalEstimate
@test rot_est.radii == radii
@test rot_est[1, 2] isa Spmt.RotationalEstimate