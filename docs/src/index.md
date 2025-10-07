```@meta
CurrentModule = SpatialMultitaper
```

# SpatialMultitaper

`SpatialMultitaper.jl` provides tools for spectral analysis of spatial point processes and random fields using multitaper methods.

## Installation
```julia
using Pkg
Pkg.add("SpatialMultitaper")
```

## Quick Start
```@example quick_start
using SpatialMultitaper, GeoStatsProcesses

import GLMakie as Mke

region = Box(Point(0, 0), Point(100, 100))
X = rand(PoissonProcess(0.01), region)
Y = rand(PoissonProcess(0.01), region)
data = spatial_data((X, Y), region)

tapers = sin_taper_family((4, 4), region)
nk = (100, 100)
kmax = (0.1, 0.1)
spec = spectra(data; tapers = tapers, nk = nk, kmax = kmax)
```

We can subset this object
```@example quick_start
spec[1,2]
```
An then plot transformations
```@example quick_start
Mke.heatmap(collect(real(spec[1,2]))...)
```

Using `collect` on any kind of estimate will return a tuple with the first D arguments being the evaluation points, and the last being the actual values. This is useful for plotting, but in that case obviously needs to be used after indexing the processes of interest.

We have other quantities that can be computed, either directly or from the spectra:
```@example quick_start
coh = coherence(spec)
Mke.heatmap(collect(abs(coh[1,2]))...)
```

```@example quick_start
rot_spec = rotational_estimate(spec)
Mke.lines(collect(real(rot_spec[1,1]))...)
```




## Contributors

```@raw html
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
```
