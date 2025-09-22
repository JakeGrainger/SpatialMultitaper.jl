#=
# Basic Estimation
=#

# This tutorial demonstrates the basic usage of the `SpatialMultitaper.jl` package for
# spectral estimation on spatial data. We will cover how to prepare your data, perform
# spectral estimation, and visualize the results.

using SpatialMultitaper, GeoStatsProcesses

import CairoMakie as Mke

# ## Preparing Your Data
# First, we need to load and prepare our spatial data. Data should be passed into
# multitapering functions as a either `PointSet` or `GeoTable` objects. If the data is
# multivariate, it should be a `Tuple` of such objects. If the data is a `GeoTable`, and the
# domain of this data is a `PointSet`, then this is interpreted as a marked point process,
# where the first column of the references information is the mark.
# If the domain is a `CartesianGrid`, then this is interpreted as a random field sampled on
# that grid. Again we take only the first column, so you should pass multivariate processes
# as separate `GeoTable` objects, even if they are recorded at the same location.

# As a simple example, say we have three independent point processes observed on the same
# region:
region = Box(Point(0, 0), Point(100, 100))
X = rand(PoissonProcess(0.01), region)
Y = rand(PoissonProcess(0.01), region)
Z = rand(PoissonProcess(0.01), region)
data = (X, Y, Z)

# We can visualise this as follows:
viz(region, color = :gray)
viz!(X, color = :black)
viz!(Y, color = :blue)
viz!(Z, color = :red)
Mke.current_figure()

# ## Estimation
# We can perform spectral estimation using the `multitaper_estimate` function. This function
# takes the `data` and a `region` on which it is oberved as inputs. In addition, we need to
# specify the tapers to use, the number of frequencies we want to compute in each dimension
# `nfreq`, and the maximum frequency in each dimension `fmax`.
tapers = sin_taper_family((4, 4), region)
spec = multitaper_estimate(data, region; tapers = tapers, nfreq = (100, 100), fmax = (1, 1))

# ## Visualising the output
# The spectral estimate is returned as a `SpectralEstimate` object. There are various
# transformations we can apply to this object. But if we want to visualise the raw output,
# we can get the frequencies and power from the fields `freq` and `power` respectively.
# The object is multidimensional, but we can index it to get the estimate between two
# processes.
spec11 = spec[1, 1]
Mke.heatmap(spec11.freq..., real.(spec11.power))
