using LinearAlgebra
using DifferentialEquations
using Plots

# Constants
M = 1.0  # Mass of the black hole
a = 0.9

observer = (20.0, deg2rad(60), 0.0, 0.0)

callback = ContinuousCallback(disc_hit_condition, disc_hit_affect!)

tspan = (0.0, 5000.0)

# Disc parameters
r_horizon = horizon(a)
r_in = isco_radius(a)       # Inner radius of the disc
r_out = 10.0                  # Outer radius of the disc

# Arrays to store trajectories
x_vals_all = []
y_vals_all = []

# Initialize array to store intersection points
disc_hits = []

r0 = observer[1]
θ0 = observer[2]
ϕ0 = observer[3]
λ0 = observer[4]

# Transforming and solving for initial conditions
x, y = 6.0, 0.0  # Example values for (x, y) in observer's image plane
x_bh, y_bh, z_bh = transform_to_bh_coords(x, y, observer, a)
r, θ, ϕ = to_boyer_lindquist(x_bh, y_bh, z_bh, a)

Σ = r^2 + a^2 * cos(θ)^2

# Initial radial, theta, phi momentum components
p_r = - (r * sqrt(r^2 + a^2) * sin(θ) * sin(θ0) * cos(ϕ - ϕ0) + (r^2 + a^2) * cos(θ) * cos(θ0))/ (Σ)
p_θ = (r * sin(θ) * cos(θ0) - sqrt(r^2 + a^2) * cos(θ) * sin(θ0) * cos(ϕ-ϕ0))/(Σ)
p_ϕ = (sin(θ0) * sin(ϕ-ϕ0))/(sqrt(r^2 + a^2) * sin(θ))

# Metric at the initial position
g = metric(r, θ)

# Solve for p^t
p_t = solve_pt(g, p_r, p_θ, p_ϕ)

# Initial momentum vector
p = [p_t, p_r, p_θ, p_ϕ]

# Update initial state vector for ODE solver
u0 = [λ0, r, θ, ϕ, p[1], p[2], p[3], p[4]]

# Calculate initial (?) Energy (E) and Angular Momentum (Lz)
E, L_z = calculate_energy_angular_momentum(g, p_t, p_ϕ)

# Calculate the impact parameter b = Lz / E
b = L_z / E

# Continue with ODE solving using the updated u0
prob = ODEProblem(intprob!, u0, tspan)
sol = solve(prob, Tsit5(), callback=callback, abstol=1e-12, reltol=1e-12, dtmax=1.0)


# Extract positions and times
t_vals = sol.t
r_vals = [sol[i][2] for i in 1:length(sol)]
θ_vals = [sol[i][3] for i in 1:length(sol)]
ϕ_vals = [sol[i][4] for i in 1:length(sol)]

# Initialize flag and variables for disc intersection
hit_disc = false
r_hit = 0.0
ϕ_hit = 0.0
p_r_hit = 0.0
p_θ_hit = 0.0
p_ϕ_hit = 0.0

# Check for intersection with the disc
for i in 1:length(θ_vals)-1
    θ1 = θ_vals[i]
    θ2 = θ_vals[i+1]

    # Check if θ crosses π/2 between θ1 and θ2
    if (θ1 - π/2) * (θ2 - π/2) <= 0
        # Interpolate to find the crossing point
        fraction = (π/2 - θ1) / (θ2 - θ1)
        r_cross = r_vals[i] + fraction * (r_vals[i+1] - r_vals[i])
        ϕ_cross = ϕ_vals[i] + fraction * (ϕ_vals[i+1] - ϕ_vals[i])

        # Normalize ϕ_cross between 0 and 2π
        ϕ_cross = mod(ϕ_cross, 2π)

        # Check if r_cross is within the disc's radial extent
        if r_in <= r_cross <= r_out
            hit_disc = true
            r_hit = r_cross
            ϕ_hit = ϕ_cross
            p_r_hit = p_r
            p_θ_hit = p_θ
            p_ϕ_hit = p_ϕ
            break  # Exit loop after first intersection
        end
    end
end

# Record the intersection if it occurs
if hit_disc
    # Store the (b, r_hit, ϕ_hit) for later use
    push!(disc_hits, (b, r_hit, ϕ_hit))
end

# Extract positions for plotting
# this is in BH frame(?) should i change to observer's?

x_vals = [sqrt(r^2 + a^2) * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
y_vals = [sqrt(r^2 + a^2) * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]

# Append trajectories
push!(x_vals_all, x_vals)
push!(y_vals_all, y_vals)

println("Number of disc hits: ", length(disc_hits))

# Use disc_hits to analyze or visualize the geodesics that hit the disc
for (b, r_hit, ϕ_hit) in disc_hits
    println("Geodesic with impact parameter b = $b hits the disc at (r, ϕ) = ($r_hit, $ϕ_hit) with p_r = $p_r_hit, p_θ = $p_θ_hit, p_ϕ = $p_ϕ_hit") 
end

#### PLOTTING ################################################################################

# Convert the disc hits to Cartesian coordinates
x_hits = [r_hit * sin(π/2) * cos(ϕ_hit) for (_, r_hit, ϕ_hit) in disc_hits]
y_hits = [r_hit * sin(π/2) * sin(ϕ_hit) for (_, r_hit, ϕ_hit) in disc_hits]

# Plot the disc boundaries
θ_values = range(0, 2π, length=500)
x_inner = [r_in * cos(θ) for θ in θ_values]
y_inner = [r_in * sin(θ) for θ in θ_values]
x_outer = [r_out * cos(θ) for θ in θ_values]
y_outer = [r_out * sin(θ) for θ in θ_values]

# Initialize the plot
pl_disc = plot(
    x_outer, y_outer,
    seriestype = :shape,
    fillcolor = :mediumorchid2,
    linecolor = :rebeccapurple,
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    title = "Disc Hits Visualization",
    label = "Disc",
    xlim = (-20, 20),
    ylim = (-20, 20),
)

# Fill the inner disc area
plot!(
    x_inner, y_inner,
    seriestype = :shape,
    fillcolor = :white,
    linecolor = :rebeccapurple,
    label = false
)

# Overlay the event horizon
circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]
plot!(
    circle_x, circle_y,
    lw = 2,
    color = :purple4,
    label = "Event horizon"
)

# Overlay the ISCO
isco_x = [r_in * cos(θ) for θ in 0:0.01:2π]
isco_y = [r_in * sin(θ) for θ in 0:0.01:2π]
plot!(
    isco_x, isco_y,
    lw = 2,
    color = :mediumslateblue,
    label = "ISCO"
)

# Overlay the hit points
scatter!(
    pl_disc,
    x_hits, y_hits,
    color = :tan1,
    markerstrokecolor = :black,
    markersize = 6,
    label = "Disc Hits",
)

# Plot the trajectory
plot!(
    pl_disc,
    x_vals, y_vals,
    lw = 1.0,
    linecolor = :lightsteelblue,
    label = false
)

# Display the plot
display(pl_disc)



###########################################################################
#   TO TRACE BACK

# r0 = r_hit
# ϕ0 = ϕ_hit
# θ0 = π/2

# tspan = (5000, 0)

# λ0 = 0.0

# p_r = -p_r_hit
# p_θ = -p_θ_hit
# p_ϕ = -p_ϕ_hit

# # Metric at the initial position
# g = metric(r0, θ0)

# # Solve for p^t
# p_t = solve_pt(g, p_r, p_θ, p_ϕ)

# norm = g[1, 1] * p_t^2 + g[2,2] * p_r^2 + g[3,3] * p_θ^2 + g[4, 4] * p_ϕ^2

# # Initial momentum vector
# p = [p_t, p_r, p_θ, p_ϕ]

# u0 = [λ0, r0, θ0, ϕ0, p[1], p[2], p[3], p[4]]

# # Continue with ODE solving using the updated u0
# prob = ODEProblem(intprob!, u0, tspan)
# sol = solve(prob, Tsit5(), callback=callback, abstol=1e-12, reltol=1e-12, dtmax=1.0)


# Reverse the momenta for time-reversed geodesic
p_r_rev = -p_r_hit
p_θ_rev = -p_θ_hit
p_ϕ_rev = -p_ϕ_hit

# Construct new initial conditions
u0_rev = [λ0, r_hit, π/2, ϕ_hit, p_t, p_r_rev, p_θ_rev, p_ϕ_rev]

# Adjust tspan for reverse time integration
tspan_rev = (0.0, -5.0)

# Solve ODE backward
prob_rev = ODEProblem(intprob!, u0_rev, tspan_rev)
sol_rev = solve(prob_rev, Tsit5(), abstol=1e-12, reltol=1e-12, dtmax=1.0)


# Extract positions and times
t_vals = sol.t
r_vals = [sol[i][2] for i in 1:length(sol)]
θ_vals = [sol[i][3] for i in 1:length(sol)]
ϕ_vals = [sol[i][4] for i in 1:length(sol)]


# Extract positions for plotting
# this is in BH frame(?) should i change to observer's?

x_vals = [sqrt(r^2 + a^2) * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
y_vals = [sqrt(r^2 + a^2) * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]

# Append trajectories
push!(x_vals_all, x_vals)
push!(y_vals_all, y_vals)

# Plot the disc boundaries
θ_values = range(0, 2π, length=500)
x_inner = [r_in * cos(θ) for θ in θ_values]
y_inner = [r_in * sin(θ) for θ in θ_values]
x_outer = [r_out * cos(θ) for θ in θ_values]
y_outer = [r_out * sin(θ) for θ in θ_values]

# Initialize the plot
pl_disc = plot(
    x_outer, y_outer,
    seriestype = :shape,
    fillcolor = :mediumorchid2,
    linecolor = :rebeccapurple,
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    title = "Disc Hits Visualization",
    label = "Disc",
    xlim = (-20, 20),
    ylim = (-20, 20),
)

# Fill the inner disc area
plot!(
    x_inner, y_inner,
    seriestype = :shape,
    fillcolor = :white,
    linecolor = :rebeccapurple,
    label = false
)

# Overlay the event horizon
circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]
plot!(
    circle_x, circle_y,
    lw = 2,
    color = :purple4,
    label = "Event horizon"
)

# Overlay the ISCO
isco_x = [r_in * cos(θ) for θ in 0:0.01:2π]
isco_y = [r_in * sin(θ) for θ in 0:0.01:2π]
plot!(
    isco_x, isco_y,
    lw = 2,
    color = :mediumslateblue,
    label = "ISCO"
)

# Overlay the hit points
scatter!(
    pl_disc,
    x_hits, y_hits,
    color = :tan1,
    markerstrokecolor = :black,
    markersize = 6,
    label = "Disc Hits",
)

# Plot the trajectory
plot!(
    pl_disc,
    x_vals, y_vals,
    lw = 1.0,
    linecolor = :lightsteelblue,
    label = false
)

# Display the plot
display(pl_disc)
