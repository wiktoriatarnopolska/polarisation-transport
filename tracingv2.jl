using LinearAlgebra
using DifferentialEquations
using Plots

# Constants
M = 1.0  # Mass of the black hole
a = 0.998

observer = (1000.0, deg2rad(60), 0.0, 0.0)

callback = ContinuousCallback(disc_hit_condition, disc_hit_affect!)

tspan = (0.0, 5000.0)

# Define a range of x values in the observer's image plane
x_values = collect(-10.0:1.0:10.0)  # x values from -10 to 10 with step size 1

# Disc parameters
r_horizon = horizon(0.9)
r_in = isco_radius(0.9)       # Inner radius of the disc
r_out = 10.0                  # Outer radius of the disc

# Arrays to store trajectories for photons that hit the disc
x_vals_hits = []
y_vals_hits = []

# Initialize array to store intersection points for photons that hit the disc
disc_hits = []

# Loop over all x values to compute photon trajectories
for x in x_values
    y = 0.0  # Could also vary y for more coverage

    # Transform to BH Cartesian coordinates
    x_bh, y_bh, z_bh = transform_to_bh_coords(x, y, observer, a)

    # Convert BH Cartesian coordinates to Boyer-Lindquist coordinates
    r, θ, ϕ = to_boyer_lindquist(x_bh, y_bh, z_bh, a)

    # Calculate initial radial, theta, and phi momentum components
    Σ = r^2 + a^2 * cos(θ)^2
    p_r = - (r * sqrt(r^2 + a^2) * sin(θ) * sin(observer[2]) * cos(ϕ - observer[3]) + (r^2 + a^2) * cos(θ) * cos(observer[2])) / Σ
    p_θ = (r * sin(θ) * cos(observer[2]) - sqrt(r^2 + a^2) * cos(θ) * sin(observer[2]) * cos(ϕ - observer[3])) / Σ
    p_ϕ = (sin(observer[2]) * sin(ϕ - observer[3])) / (sqrt(r^2 + a^2) * sin(θ))

    # Metric at the initial position
    g = metric(r, θ)

    # Solve for p^t using the null condition
    p_t = solve_pt(g, p_r, p_θ, p_ϕ)

    # Initial momentum vector
    p = [p_t, p_r, p_θ, p_ϕ]

    # Update initial state vector for ODE solver
    u0 = [observer[4], r, θ, ϕ, p[1], p[2], p[3], p[4]]

    # Define the ODE problem for the photon
    prob = ODEProblem(intprob!, u0, tspan)

    # Solve the ODE
    sol = solve(prob, Tsit5(), callback=callback, abstol=1e-12, reltol=1e-12, dtmax=1.0)

    # Extract radial and azimuthal values for plotting
    r_vals = [sol[i][2] for i in 1:length(sol)]
    θ_vals = [sol[i][3] for i in 1:length(sol)]
    ϕ_vals = [sol[i][4] for i in 1:length(sol)]

    # Initialize flag and variables for disc intersection
    hit_disc = false
    r_hit = 0.0
    ϕ_hit = 0.0

    # Check for intersection with the disc
    for i in 1:length(θ_vals) - 1
        θ1 = θ_vals[i]
        θ2 = θ_vals[i + 1]

        # Check if θ crosses π/2 between θ1 and θ2
        if (θ1 - π/2) * (θ2 - π/2) <= 0
            # Interpolate to find the crossing point
            fraction = (π/2 - θ1) / (θ2 - θ1)
            r_cross = r_vals[i] + fraction * (r_vals[i + 1] - r_vals[i])
            ϕ_cross = ϕ_vals[i] + fraction * (ϕ_vals[i + 1] - ϕ_vals[i])
            
            # Normalize ϕ_cross between 0 and 2π
            ϕ_cross = mod(ϕ_cross, 2π)

            # Check if r_cross is within the disc's radial extent
            if r_in <= r_cross <= r_out
                hit_disc = true
                r_hit = r_cross
                ϕ_hit = ϕ_cross
                break  # Exit loop after first intersection
            end
        end
    end

    # Record the intersection and store the trajectory if it occurs
    if hit_disc
        # Store the (x, r_hit, ϕ_hit) for later use
        push!(disc_hits, (x, r_hit, ϕ_hit))

        # Extract positions for plotting
        x_vals = [sqrt(r^2 + a^2) * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
        y_vals = [sqrt(r^2 + a^2) * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]


        # Append trajectories for photons that hit the disc
        push!(x_vals_hits, x_vals)
        push!(y_vals_hits, y_vals)
    end
end

# Print the number of disc hits
println("Number of disc hits: ", length(disc_hits))

# Use disc_hits to analyze or visualize the geodesics that hit the disc
for (b, r_hit, ϕ_hit) in disc_hits
    println("Geodesic with impact parameter b = $b hits the disc at (r, ϕ) = ($r_hit, $ϕ_hit)")
end

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
    legend=false
)

# Fill the inner disc area
plot!(
    x_inner, y_inner,
    seriestype = :shape,
    fillcolor = :white,
    linecolor = :rebeccapurple,
    label = false
)

# Overlay the hit points
for (x, r_hit, ϕ_hit) in disc_hits
    x_hit = r_hit * cos(ϕ_hit)
    y_hit = r_hit * sin(ϕ_hit)
    scatter!(
        pl_disc,
        [x_hit], [y_hit],
        color = :tan1,
        markerstrokecolor = :black,
        markersize = 6,
        label = "Disc Hit (x=$x)"
    )
end

# Plot the geodesic paths of only those photons that hit the disc
for i in 1:length(x_vals_hits)
    plot!(
        pl_disc,
        x_vals_hits[i], y_vals_hits[i],
        lw = 1.0,
        linecolor = :lightsteelblue,
        label = false
    )
end

# Overlay the event horizon
circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]
plot!(
    pl_disc,
    circle_x, circle_y,
    lw = 2,
    color = :purple4,
    label = "Event horizon"
)

# Overlay the ISCO
isco_x = [r_in * cos(θ) for θ in 0:0.01:2π]
isco_y = [r_in * sin(θ) for θ in 0:0.01:2π]
plot!(
    pl_disc,
    isco_x, isco_y,
    lw = 2,
    color = :mediumslateblue,
    label = "ISCO"
)

# Display the plot
display(pl_disc)
