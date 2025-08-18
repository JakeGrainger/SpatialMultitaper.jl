region = Box(Point(0, 0), Point(3, 3))
pattern = PointSet([Point(0, 0), Point(1, 1)])
pattern2 = PointSet([Point(0.3, 0.2), Point(0.8, 0.4), Point(0.5, 0.5)])
pattern3 = PointSet([Point(0.1, 0.1), Point(0.2, 0.3), Point(0.4, 0.6)])
# griddata = georef((rf = rand(6 * 6),), CartesianGrid((0, 0), (3, 3), dims = (6, 6)))
tapers = sin_taper_family((3, 3), region)
nfreq = (10, 10)
fmax = (2, 2)
data = (pattern, pattern2, pattern3)
radii = 0.3:0.1:1.0


results = C_function(pattern, region, radii; tapers = tapers, nfreq = nfreq, fmax = fmax)
@test results isa Spmt.CFunction
results = K_function(pattern, region, radii; tapers = tapers, nfreq = nfreq, fmax = fmax)
@test results isa Spmt.KFunction
# results = L_function(pattern, region, radii; tapers = tapers, nfreq = nfreq, fmax = fmax)
# @test results isa Spmt.LFunction

results = C_function(data, region, radii; tapers = tapers, nfreq = nfreq, fmax = fmax)
@test results isa Spmt.CFunction
@test results[1, 2].C_function isa Vector{Float64}
@test results[1, 2].C_function == Spmt.getestimate(results)[(1, 2)]

results = K_function(data, region, radii; tapers = tapers, nfreq = nfreq, fmax = fmax)
@test results isa Spmt.KFunction
# results = L_function(data, region, radii; tapers = tapers, nfreq = nfreq, fmax = fmax)
# @test results isa Spmt.LFunction

results =
    partial_C_function(data, region, radii; tapers = tapers, nfreq = nfreq, fmax = fmax)
@test results isa Spmt.PartialCFunction

results =
    partial_K_function(data, region, radii; tapers = tapers, nfreq = nfreq, fmax = fmax)
@test results isa Spmt.PartialKFunction

# results =
#     partial_L_function(data, region, radii; tapers = tapers, nfreq = nfreq, fmax = fmax)
# @test results isa Spmt.PartialLFunction

results = paircorrelation_function(
    data,
    region,
    radii;
    tapers = tapers,
    nfreq = nfreq,
    fmax = fmax,
)
@test results isa Spmt.PairCorrelationFunction

results = paircorrelation_function_direct(
    data,
    region,
    radii;
    tapers = tapers,
    nfreq = nfreq,
    fmax = fmax,
)
@test results isa Spmt.PairCorrelationFunction

results = partial_paircorrelation_function(
    data,
    region,
    radii;
    tapers = tapers,
    nfreq = nfreq,
    fmax = fmax,
)
@test results isa Spmt.PartialPairCorrelationFunction

results = partial_paircorrelation_function_direct(
    data,
    region,
    radii;
    tapers = tapers,
    nfreq = nfreq,
    fmax = fmax,
)
@test results isa Spmt.PartialPairCorrelationFunction
