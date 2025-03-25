mt_est = make_simple_example()
x = partial_complex_coherence(mt_est)
@test x isa Spmt.PartialComplexCoherence
@test x.freq == mt_est.freq
