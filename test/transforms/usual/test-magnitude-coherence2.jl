mt_est = make_simple_example()
coh2 = magnitude_coherence2(mt_est)
@test coh2 isa Spmt.MagnitudeCoherence2
@test coh2.freq == mt_est.freq
coh = magnitude_coherence(mt_est)
@test coh2[1,2].coherence ≈ coh[1,2].coherence.^2