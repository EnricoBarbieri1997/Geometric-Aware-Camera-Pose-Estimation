using CylindersBasedCameraResectioning.Geometry: get_tangentpoints_circle_point, Circle
using CairoMakie

# Example data
center = [0.8, 1.2, 0.0]
radius = 2.0
external_point = [4.0, 3.0, 0.0]

circle = Circle(center, radius, [0, 0, 1])

# Get tangent points
tp1, tp2 = get_tangentpoints_circle_point(circle, external_point)

# Plot
fig = Figure()

ax = Axis(fig[1, 1], aspect = DataAspect())

# Circle
circle = Makie.Circle(Makie.Point2f0(center[1:2]...), radius)
poly!(ax, circle, color = :white, strokewidth = 2, strokecolor = :black)

# External point
scatter!(ax, [external_point[1]], [external_point[2]], color = :red, markersize = 10)

# Tangent points
scatter!(ax, [tp1[1], tp2[1]], [tp1[2], tp2[2]], color = :blue, markersize = 8)

# Tangent lines
lines!(ax, [external_point[1], tp1[1]], [external_point[2], tp1[2]], color = :green, linewidth = 2)
lines!(ax, [external_point[1], tp2[1]], [external_point[2], tp2[2]], color = :green, linewidth = 2)

fig