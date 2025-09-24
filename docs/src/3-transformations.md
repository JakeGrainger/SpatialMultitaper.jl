# Transformations

There are various transformations that one can compute from estimates of the spectral density function. Each of them can either be applied to a `Matrix`, or a `SpectralEstimate`.
As a running example, consider a simulated example of three independent Poisson processes:
```@example poisson_example
import GLMakie as Mke
using SpatialMultitaper, GeoStatsProcesses
region = Box(Point(0,0), Point(100, 100))
X = rand(PoissonProcess(0.01), region)
Y = rand(PoissonProcess(0.01), region)
Z = rand(PoissonProcess(0.01), region)
data = (X, Y, Z)
fig = Mke.Figure()
ax = Mke.Axis(fig[1,1], aspect = 1)
viz!(ax, X, label = "X", color = :blue, pointmarker = 'X')
viz!(ax, Y, label = "Y", color = :orange, pointmarker = 'Y')
viz!(ax, Z, label = "Z", color = :purple, pointmarker = 'Z')
fig
```

We can then estimate the spectral density function
```@example poisson_example
tapers = sin_taper_family((4,4), region)
sdf = multitaper_estimate(data, region; tapers = tapers, nfreq = (100,100), fmax = (1,1))
```

## Coherence

Coherence is a frequency domain notion of correlation.
This quantity is complex valued, and so there are a variety of ways to think about it.
In particular, one often considers some combination of complex coherence, magnitude coherence, magnitude square coherence and group delay.
Complex coherence is defined as
```math
    c_{X,Y}(k) = \frac{S_{XY}(k)}{\sqrt{S_{XX}(k)S_{YY}(k)}}.
```
Coherence is then the magnitude of this quantity, and magnitude square coherence is the square of this quantity, and group delay is its complex argument.

## Phase

## Partial coherence

## Partial phase

## Partial spectral density
