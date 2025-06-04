mt_est = make_simple_example()
x = partial_phase(mt_est)
@test x isa Spmt.PartialPhase
@test x.freq == mt_est.freq
