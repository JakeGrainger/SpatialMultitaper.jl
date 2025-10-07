```@meta
EditURL = "../../literate/tutorials/k_function.jl"
```

# K Function

````@example k_function
using SpatialMultitaper, GeoStatsProcesses

import GLMakie as Mke
````

The K function is a second-order summary statistic for point processes. It is defined as
the expected number of further points within a distance `r` of an arbitrary point, divided
by the intensity of the process. For a homogeneous Poisson process, this is simply `πr^2`
in two dimensions. The K function can be estimated using the `k_function` function.

````@example k_function
region = Box(Point(0, 0), Point(100, 100))
X = rand(PoissonProcess(0.01), region)
data = spatial_data(X, region)
tapers = sin_taper_family((4, 4), region)
nk = (100, 100)
kmax = (0.1, 0.1)
radii = 0:0.5:30
kfun = k_function(data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)
````

We can then visualise the K function using Makie

````@example k_function
Mke.lines(kfun)
````

## L functions
The L function is a transformation of the K function, defined as
`L(r) = sqrt(K(r)/π)`. This is usually easier to visualise, as for a homogeneous Poisson
process, `L(r) = r`. We can compute the L function using the `l_function` function.

````@example k_function
lfun = l_function(kfun)
````

and again plot this using Makie

````@example k_function
Mke.lines(lfun)
````

we can also plot the difference `L(r) - r`, which is often used to assess clustering or
inhibition in the point process.

````@example k_function
Mke.lines(centered_l_function(lfun))
````

these functions can be computed from each other or directly from the data

````@example k_function
lfun2 = l_function(data, radii = radii, nk = nk, kmax = kmax, tapers = tapers)
lfun.value == lfun2.value
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

