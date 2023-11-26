
f = Figure(size = (1000, 700))
ax = Axis3(f[1, 1], title = "Cylinders")

colors = [:red, :green, :blue, :yellow]

for i in 1:numberOfCylinders
    cyl1 = CylinderMesh(Point(0.0,0.0,0.0), Point(0.0,0.0,1.0), 1.0)
    scale!(cyl1, radiuses[i][1], radiuses[i][1])
    # cyl1 = cyl1 |> Stretch(radiuses[i][1], radiuses[i][1], 1)
    mesh!(ax, cyl1, color = colors[i], model = transforms[i]) #
end

f