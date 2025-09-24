# Diagrams

This diagram shows how the different statistical measures and spatial functions in SpatialMultitaper.jl relate to each other. Each estimate type can be computed in both **Full** and **Partial** versions using the same transformation functions.

```mermaid
graph TD
    S[Spectra] --> Coh[Coherence]
    S --> P[Phase]
    Coh --> M[Magnitude coherence]
    Coh --> M2[Magnitude squared coherence]
    Coh --> P
    S --> C[C function]
    C --> K[K function]
    K --> L[L function]
    L --> cL[centered L function]
```

**Note**: All transformation relationships shown above work identically for both Full and Partial estimate types.