integral, err = quadgk(k -> Spmt.sphere_weight(1.4, k, Val{2}()), 0.5, 0.7)
@test integral ≈ Spmt.iso_weight(1.4, 0.6, 0.2, Val{2}())
