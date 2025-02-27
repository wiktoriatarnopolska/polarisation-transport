using Plots

N_r = 50
N_phi = 50

rin = isco_radius(0.9)
rout = 10

# Geometric radial grid: Logarithmically spaced
r_grid_log = range(log10(rin), stop=log10(rout), length=N_r)
r_grid = 10 .^ r_grid_log  # Convert back to linear space

# Option: Linear radial grid (risk of losing signal from the ISCO)
#r_grid = range(rin, rout, length=N_r)  # Radial grid in units of r_g

phi_grid = range(0, 2π, length=N_phi)  # Azimuthal grid from 0 to 2π

# Convert polar coordinates (r, φ) to Cartesian coordinates (x, y) for plotting
x_vals, y_vals = Float64[], Float64[]

# cross-check with references
# Impact parameters
α_vals, β_vals = Float64[], Float64[]

for r in r_grid
    for φ in phi_grid
        x = r * cos(φ)
        y = r * sin(φ)
        push!(x_vals, x)
        push!(y_vals, y)

        # Impact parameters
        α = - r * sin(φ)
        β = r * cos(φ) * sin(θobs)
        push!(α_vals, α)
        push!(β_vals, β)
    end
end


# Plot the grid points
scatter(x_vals, y_vals, xlabel="x (r_g)", ylabel="y (r_g)", 
        title="Sampled Disk Grid", aspect_ratio=:equal, label = "x and y coordinates")

# Plot impact parameters
scatter!(α_vals, β_vals, 
        #xlabel = "α (r_g)", ylabel = "β (r_g)",
        #title="Impact parameters", 
        label = "impact parameters α and β",
        aspect_ratio=:equal)
# savefig("grid.png")
