mt_est = make_simple_example()
coh = complex_coherence(mt_est)
@test coh isa Spmt.ComplexCoherence
@test coh.freq == mt_est.freq
@test coh[1,1].coherence[1] ≈ complex(1.0)
@test coh[1,2].coherence[1] ≈ mt_est[1,2].power[1] / sqrt(mt_est[1,1].power[1] * mt_est[2,2].power[1])