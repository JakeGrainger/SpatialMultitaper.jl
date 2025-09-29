# Multitapering

Multitapering is a popular technique for estimating the spectral density function of a signal in time. Multitapering can also be extended to spatial processes, including both random fields [hanssen1997multidimensional](@citep), point processes [rajala2023what](@citep) and multivariate processes which are a mixture of the two [grainger2025spectral](@citep).

## Data format

Data should be passed into multitapering functions as a either `PointSet` or `GeoTable` objects.
If the data is multivariate, it should be a `Tuple` of such objects.
If the data is a `GeoTable`, and the domain of this data is a `PointSet`, then this is interpreted as a marked point process, where the first column of the references information is the mark. If the domain is a `CartesianGrid`, then this is interpreted as a random field sampled on that grid. Again we take only the first column, so you should pass multivariate processes as separate `GeoTable` objects, even if they are recorded at the same location.

## Multivariate processes

```@docs
spectra
```

## Mathematical background

# References

```@bibliography
Pages = ["2-multitapering.md"]
```