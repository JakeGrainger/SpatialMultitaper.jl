mt_est = make_simple_example()
coh2 = partial_magnitude_coherence2(mt_est)
@test coh2 isa Spmt.PartialMagnitudeCoherence2
@test coh2.freq == mt_est.freq
coh = partial_magnitude_coherence(mt_est)
@test coh2[1,2].partial_coherence ≈ coh[1,2].partial_coherence.^2