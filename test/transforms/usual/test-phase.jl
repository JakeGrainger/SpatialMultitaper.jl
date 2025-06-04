mt_est = make_simple_example()
x = phase(mt_est)
@test x isa Spmt.Phase
@test x.freq == mt_est.freq
