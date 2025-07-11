using SpatialMultitaper
import CairoMakie as Mke

c = 0.3
region = PolyArea([
    (0, 1 - c),
    (c, 1),
    (1 - c, 1),
    (1, 1 - c),
    (1, c),
    (1 - c, 0),
    (c, 0),
    (0, c),
])
viz(region)

tapers, λ = make_tapers(
    region,
    bandwidth = 3.5,
    space_res = 400,
    freq_res = 500,
    freq_generate_res = 1000,
)

pattern = georef((marks = [1, 2.4],), PointSet([Point(0, 0), Point(1, 1)]))
griddata = georef((rf = rand(200 * 200),), CartesianGrid((0, 0), (1, 1), dims = (200, 200)))
griddata2 = georef((rf = rand(10 * 10),), CartesianGrid((0, 0), (1, 1), dims = (10, 10)))

data = (pattern, griddata, griddata2)
tapers_on_grids =
    SpatialMultitaper.tapers_on_grid.(
        Ref(tapers),
        SpatialMultitaper.get_grid.(data),
        freq_res = 500,
    );


normalisation, conc =
    SpatialMultitaper.taper_checks(data, tapers, 3.5, freq_res = 500, conc_res = 200)
normalisation
[conc[i, i, :, :] for i in axes(conc, 1)][1]

##
n_space = 400
Mke.lines(range(-10, 10, 100), x -> sinc(x * (1 / n_space))^4)
Mke.ylims!(Mke.current_axis(), 0, 1)
Mke.current_figure()


##
fig = Mke.Figure()
ax = [Mke.Axis(fig[i, j], aspect = 1) for i = 1:3, j = 1:4]
for i in eachindex(tapers)
    Mke.heatmap!(
        ax[i],
        range(-5, 5, 100),
        range(-5, 5, 100),
        (x, y) -> abs2(taper_ft(tapers_on_grids[3][i], x, y)),
        colorrange = (0, 0.2),
    )
    viz!(ax[i], boundary(Ball(Point(0, 0), 3.5)), color = :black)
end
fig

##
fig = Mke.Figure()
ax = [Mke.Axis(fig[1, i], aspect = 1) for i = 1:3]
for i = 1:3
    Mke.heatmap!(
        ax[i],
        range(-5, 5, 100),
        range(-5, 5, 100),
        (x, y) -> abs2(taper_ft(tapers_on_grids[i][1], x, y)),
        colorrange = (0, 0.3),
    )
    viz!(ax[i], boundary(Ball(Point(0, 0), 3.5)), color = :black)
end
fig

check(x, y) =
    abs2(taper_ft(tapers_on_grids[1][1], x, y)) /
    abs2(taper_ft(tapers_on_grids[3][1], x, y))

##
fig = Mke.Figure()
ax = Mke.Axis(fig[1, 1])
for i = 1:3
    Mke.lines!(ax, 0:0.01:3.5, (x) -> abs(taper_ft(tapers_on_grids[i][1], x, 0.0)))
end
fig

##
