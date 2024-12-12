mt_est = make_simple_example()
x = partial_magnitude_coherence(mt_est)
@test x isa Spmt.PartialMagnitudeCoherence
@test x.freq == mt_est.freq
