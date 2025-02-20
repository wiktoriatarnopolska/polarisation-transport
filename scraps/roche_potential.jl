using Plots

# Define the Roche potential function
function roche_potential(x, y, q)
    r1 = sqrt((x + q)^2 + y^2)  # Distance to the primary star
    r2 = sqrt((x - 1 + q)^2 + y^2)  # Distance to the secondary star
    return -(1-q)/r1 - q/r2 - 0.5*(x^2 + y^2)
end

M1 = 10.0      # Mass of BH in solar masses
M2 = 2.5       # Mass of companion
q = M2 / M1    # Mass ratio

a = 1

x_range = (-2:0.05:2)*a
y_range = (-2:0.05:2)*a 

x = collect(x_range)
y = collect(y_range)
z = [roche_potential(xi, yi, q) for yi in y, xi in x]

x_grid = x[1:3:end]
y_grid = y[1:3:end]
z_grid = [roche_potential(xi, yi, q) for yi in y_grid, xi in x_grid]

wireframe(x_grid, y_grid, z_grid)

plot!(x, y, z, st=:surface, color=:linear_worb_100_25_c53_n256,
xlabel=L"x \, (R_\textrm{⊙})", zlabel=L"\textrm{Potential} \, (\textrm{dimensionless})",
camera=(38, 50), zlims = (-6, 0), fillalpha = 0.9, legend=:left, frame=:box)  # Azimuth and Elevation angles
annotate!(0, 5, -16, text(L"y \, (R_\textrm{⊙})", 12, :black))
