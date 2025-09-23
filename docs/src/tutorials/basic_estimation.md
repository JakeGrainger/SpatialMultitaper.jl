```@meta
EditURL = "../../literate/tutorials/basic_estimation.jl"
```

# Basic Estimation

Multitapering is a popular technique for estimating the spectral density function of a
signal in time[thomson1982spectrum](@citep). Multitapering can also be extended to spatial processes, including both
random fields [hanssen1997multidimensional](@citep),
point processes [rajala2023what](@citep) and multivariate processes which are a mixture of
the two [grainger2025spectral](@citep).

This tutorial demonstrates the basic usage of the `SpatialMultitaper.jl` package for
spectral estimation on spatial data. We will cover how to prepare your data, perform
spectral estimation, and visualize the results.

````@example basic_estimation
using SpatialMultitaper, GeoStatsProcesses

import CairoMakie as Mke
````

Set default aspect ratio for plots

````@example basic_estimation
Mke.set_theme!(; Axis = (aspect = Mke.DataAspect(),))
````

## Preparing Your Data
First, we need to load and prepare our spatial data. Data should be passed into
multitapering functions as a either `PointSet` or `GeoTable` objects. If the data is
multivariate, it should be a `Tuple` of such objects. If the data is a `GeoTable`, and the
domain of this data is a `PointSet`, then this is interpreted as a marked point process,
where the first column of the references information is the mark.
If the domain is a `CartesianGrid`, then this is interpreted as a random field sampled on
that grid. Again we take only the first column, so you should pass multivariate processes
as separate `GeoTable` objects, even if they are recorded at the same location.

As a simple example, say we have three independent point processes observed on the same
region:

````@example basic_estimation
region = Box(Point(0, 0), Point(100, 100))
X = rand(PoissonProcess(0.01), region)
Y = rand(PoissonProcess(0.01), region)
data = (X, Y)
````

We can visualise this as follows:

````@example basic_estimation
viz(boundary(region), color = :gray)
viz!(X, color = :black)
viz!(Y, color = :red)
Mke.current_figure()
````

## Estimation
We can perform spectral estimation using the `multitaper_estimate` function. This function
takes the `data` and a `region` on which it is oberved as inputs. In addition, we need to
specify the tapers to use, the number of frequencies we want to compute in each dimension
`nfreq`, and the maximum frequency in each dimension `fmax`.

````@example basic_estimation
tapers = sin_taper_family((4, 4), region)
nfreq = (100, 100)
fmax = (0.1, 0.1)
spec = multitaper_estimate(data, region; tapers = tapers, nfreq = nfreq, fmax = fmax)
````

## Visualising the output
The spectral estimate is returned as a `SpectralEstimate` object. There are various
transformations we can apply to this object. But if we want to visualise the raw output,
we can get the frequencies and power from the fields `freq` and `power` respectively.
The object is multidimensional, but we can index it to get the estimate between two
processes.

````@example basic_estimation
spec11 = spec[1, 1]
Mke.heatmap(spec11.freq..., real.(spec11.power))
````

## A more interesting example
Let's consider a more interesting example, where we have a bivariate process with some
cross-correlation. We can simulate a very simple process by

````@example basic_estimation
region = Box(Point(0, 0), Point(100, 100))
shift = 10.0
big_region = Box(Point(-shift, -shift), Point(100 + shift, 100 + shift))
X = rand(PoissonProcess(0.01), big_region)
Y = Translate(shift, shift)(X)
data = mask.((X, Y), Ref(region))

tapers = sin_taper_family((4, 4), region)
nfreq = (100, 100)
fmax = (0.1, 0.1)
spec = multitaper_estimate(data, region; tapers = tapers, nfreq = nfreq, fmax = fmax)

Mke.heatmap(collect(real(spec[1, 2]))...)

example_coh = magnitude_coherence(spec)
example_phase = phase(spec)
fig = Mke.Figure()
Mke.heatmap(fig[1, 1], collect(example_coh[1, 2])..., colorrange = (0, 1))
Mke.heatmap(fig[1, 2], collect(example_phase[1, 2])..., colorrange = (-π, π))
fig
````

## References

```@bibliography
Pages = ["basic_estimation.md"]
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

