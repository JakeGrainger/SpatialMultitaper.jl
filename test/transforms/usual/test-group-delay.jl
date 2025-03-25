mt_est = make_simple_example()
x = group_delay(mt_est)
@test x isa Spmt.GroupDelay
@test x.freq == mt_est.freq
