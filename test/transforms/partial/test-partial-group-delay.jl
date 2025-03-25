mt_est = make_simple_example()
x = partial_group_delay(mt_est)
@test x isa Spmt.PartialGroupDelay
@test x.freq == mt_est.freq
