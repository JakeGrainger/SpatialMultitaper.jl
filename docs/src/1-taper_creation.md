# Tapers

Tapers are an important part of multitapering. They are used both to reduce bias from edge effects, and to reduce variance by providing "pseudo-replications".
In `SpatialMultitaper`, there are a number of dedicated structs to represent the notion of a taper, as well as many convnience functions to create them.
It is important to distinguish between tapers for continuous processes, and tapers for discrete processes. However, we shall begin the description for the continuous case, as one can often subsample these to obtain legitimate discrete tapers. See the [mathematical background](#mathematical-background) and [developer notes](#developer-notes) for more information.

## Creating tapers on a box

If the region of interest is a box, then one can easily extend a family of univariate tapers by taking the family of their cross products.
One dimensional families on the unit interval are easily extended to a general interval. In particular, say that we have a taper ``\tilde h:[0,1]\rightarrow \mathbb{R}``, then we may define ``h:[a,b]\rightarrow \mathbb{R}`` by ``h(x) = \tilde h((x-a)/(b-a)) / \sqrt{b-a}``.
The simplest family of continuous tapers on a box are the sin tapers [riedel1995minimum](@citep).
```@docs 
sin_taper_family
```

## Creating tapers on arbitrary regions

```@example octogon_taper
using SpatialMultitaper # hide
import CairoMakie as Mke # hide
```

We can also construct families of tapers for other regions of interest, using the methodology proposed by [simons2011spatiospectral](@cite).
For example, say we are interested in an octagonal region

```@example octogon_taper
c = 0.3
region = PolyArea([
  (0,1-c), (c,1), (1-c,1), (1,1-c), (1,c), (1-c,0), (c,0), (0,c)
])
viz(region)
```

We can construct discrete tapers for this region using the `make_taper` function:
```@example octogon_taper
    tapers, concentrations = make_tapers(region, bandwidth = 3.5)
```
In this case, we obtain 12 well concentrated tapers within this bandwidth. We can visualize the tapers:
```@example octogon_taper
    fig = Mke.Figure()
    ax = [Mke.Axis(fig[i,j], aspect = 1) for i in 1:3, j in 1:4]
    for i in eachindex(tapers)
        Mke.heatmap!(ax[i], range(0,1,100), range(0,1,100), (x,y) -> tapers[i](x,y), colorrange = (-3,3))
        viz!(ax[i], boundary(region), color = :black)
    end
    Mke.Colorbar(fig[1:3,5], colorrange = (-3,3))
    fig
```

Note that this process can be slow for complex regions, or if high resolution grids are used.

```@docs
make_tapers
```

## Exploring the tapers

Given a constructed family of tapers, one can index into them to obtain one single taper.
```@example octogon_taper
    taper = tapers[1]
```

We can evaluate the taper at a given point by calling it as a function:
```@example octogon_taper
    taper(0.5, 0.5)
```
Alternatively, we can pass a point or a tuple
```@example octogon_taper
    taper(Point(0.5, 0.5))
```
```@example octogon_taper
    taper((0.5, 0.5))
```
We can also consider the Fourier transform of the taper, which is a function of wavenumber.
In general, we would like the Fourier transform to be well concentrated around the origin, and to decay quickly away from it.
When constructing the general tapers, we specified a bandwidth, within which we expect the Fourier transform to be well concentrated (in a ball around zero with a radius of the bandwidth).
We can evaluate the Fourier transform of the taper at a given wavenumber with the `taper_ft` function:
```@example octogon_taper
    taper_ft(taper, 0, 0)
```
So plotting the `abs2` of Fourier transforms to check the concentrations we see that
```@example octogon_taper
    fig = Mke.Figure()
    ax = [Mke.Axis(fig[i,j], aspect = 1) for i in 1:3, j in 1:4]
    for i in eachindex(tapers)
        Mke.heatmap!(ax[i], range(-5,5,100), range(-5,5,100), (x,y) -> abs2(taper_ft(tapers[i],x,y)))
        viz!(ax[i], boundary(Ball(Point(0,0), 3.5)), color = :black)
    end
    fig
```

## Mathematical background

## Developer notes

## References

```@bibliography
Pages = ["1-taper_creation.md"]
```