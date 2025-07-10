for vector_of_processes in [false, true]
    mt_est = make_simple_example(vector_of_processes = vector_of_processes)
    coh = magnitude_coherence(mt_est)
    @test coh isa Spmt.MagnitudeCoherence
    @test coh.freq == mt_est.freq
    coh_complex = complex_coherence(mt_est)
    @test map(x -> abs.(x), coh_complex.coherence) == coh.coherence
end