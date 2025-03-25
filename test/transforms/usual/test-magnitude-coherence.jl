mt_est = make_simple_example()
coh = magnitude_coherence(mt_est)
@test coh isa Spmt.MagnitudeCoherence
@test coh.freq == mt_est.freq
coh_complex = complex_coherence(mt_est)
@test map(x->abs.(x), coh_complex.coherence) == coh.coherence