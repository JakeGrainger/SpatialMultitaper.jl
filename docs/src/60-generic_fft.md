# Generic ffts

As part of the library, we allow the user to request there own output wavenumbers.
They can do this by specifying the highest wavenumber they want, and the number of wavenumbers.
The wavenumbers will then be given back in the format used by `FFTW`, but shifted. In particular,
```@example
using FFTW
n = 10
highest_freq = 2.5
freq = fftshift(fftfreq(n, 2highest_freq))
```
Let `n` be the number of frequencies and `N` be the highest frequency. 
We construct frequencies of the form
```math
    \{-\lceil n/2\rceil, \ldots, \lfloor n/2\rfloor - 1\} \times \frac{2N}{n}
```

This can be recovered by using
```@example
    using SpatialMultitaper
    n = 10
    N = 2.5
    SpatialMultitaper._choose_frequencies_1d(n, N)
```

We utilise both standard Fast Fourier Transforms (FFTs) and Non-Uniform Fast Fourier Transforms (NUFFTs)

## NUFFTs
The standard NUFFTs as implemented by `FINUFFT` which we need are the first type, which from irregular inputs to a regular output.
If the user request `ms` outputs to be given, then they get back the Fourier transform at `z/2pi` for all integers `z` between `-ms/2` and `(ms-1)/2` inclusive.
The assumption is that the original data is contained in `[-3pi, 3pi)`.

Currently, we transform data to be on `[-pi, pi)`.
Multiplying by a scaling in one domain always divides in the other domain.
Consider data recorded in the interval `[a,b)`.
Say that one such point is denoted by `x`.
We then construct a new point, say `y`, by computing
```math
    c = (b-a)/2
    y = 2\pi (x - c) / l
```
for some choice of scaling `l`.
Now the output of an NUFFT (assuming `y` is now in `[-3pi, 3pi)`) is `lz` for all integers `z` between `-ms/2` and `(ms-1)/2` inclusive.
So given a desired collection of wavenumbers, we need to pick `l` and `ms` so that we could recover the wavenumber we are interested in from the output wavenumber, and so that all points 'y' are contained in `[-pi, pi)`.
