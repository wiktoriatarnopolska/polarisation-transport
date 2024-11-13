using LinearAlgebra
using DifferentialEquations
using Plots

# Constants
M = 1.0  # Mass of the black hole
a = 0.9

observer = (1000.0, deg2rad(60), 0.0, 0.0)

callback = ContinuousCallback(disc_hit_condition, disc_hit_affect!)

tspan = (0.0, 5000.0)

# Define a range of impact parameters
b_values = collect(-10.0:1:10.0)

# Disc parameters
r_horizon = horizon(0.9)
r_in = isco_radius(0.9)       # Inner radius of the disc
r_out = 10.0                  # Outer radius of the disc

# Define the ODE function using an in-place form
function intprob!(du, u, p, λ)
    # Current position and velocity
    x = u[1:4]
    v = u[5:8]

    r_val, θ_val = x[2], x[3]

    # Initialize du
    du .= 0.0

    # Terminate the integration if r crosses the event horizon
    if r_val <= r_horizon
        du[1:4] .= v
        return
    end

    # Compute analytical Christoffel symbols
    Γ = compute_christoffel_analytical(r_val, θ_val)

    # Compute the derivatives of velocity
    dv = zeros(4)
    for μ in 1:4
        sum_ = 0.0
        for ν in 1:4
            for λ in 1:4
                sum_ += Γ[μ, ν, λ] * v[ν] * v[λ]
            end
        end
        dv[μ] = -sum_
    end

    # Assign derivatives to du
    du[1:4] .= v
    du[5:8] .= dv
end

# Arrays to store trajectories
x_vals_all = []
y_vals_all = []

# Initialize array to store intersection points
disc_hits = []

for b in b_values
    # Compute initial position in Cartesian coordinates
    r0 = observer[1]
    θ0 = observer[2]
    ϕ0 = observer[3]
    λ0 = observer[4]

    # Initial velocities in spherical coordinates
    v_r = -1.0                 # Photon moving inward
    v_θ = 0.0                  # Equatorial plane
    v_ϕ = b / (r0^2)           # Angular velocity from impact parameter.

    # Compute the metric at the initial position
    g0 = metric(r0, θ0)
    g_tt = g0[1,1]
    g_tϕ = g0[1,4]
    g_rr = g0[2,2]
    g_θθ = g0[3,3]
    g_ϕϕ = g0[4,4]

    # Solve quadratic equation for v_t
    A = g_tt
    B = g_tϕ * v_ϕ
    C = g_rr * v_r^2 + g_θθ * v_θ^2 + g_ϕϕ * v_ϕ^2

    Δ = (2 * B)^2 - 4 * A * C
    if Δ < 0
        error("No real solution for v_t. Check initial conditions.")
    end

    v_t = (-2 * B - sqrt(Δ)) / (2 * A)  # Negative root for future-directed motion

    # Complete the initial velocity vector
    v = [v_t; v_r; v_θ; v_ϕ]

    # Initial state vector
    u0 = [λ0, r0, θ0, ϕ0, v[1], v[2], v[3], v[4]]

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
    x_vals = [r * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
    y_vals = [r * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]

    # Append trajectories
    push!(x_vals_all, x_vals)
    push!(y_vals_all, y_vals)
end

println("Number of disc hits: ", length(disc_hits))

# Now you can use disc_hits to analyze or visualize the geodesics that hit the disc
for (b, r_hit, ϕ_hit) in disc_hits
    println("Geodesic with impact parameter b = $b hits the disc at (r, ϕ) = ($r_hit, $ϕ_hit)")
end

# #### PLOTTING ################################################################################

# # Convert the disc hits to Cartesian coordinates
# x_hits = [r_hit * sin(π/2) * cos(ϕ_hit) for (_, r_hit, ϕ_hit) in disc_hits]
# y_hits = [r_hit * sin(π/2) * sin(ϕ_hit) for (_, r_hit, ϕ_hit) in disc_hits]

# # Plot the disc boundaries
# θ_values = range(0, 2π, length=500)
# x_inner = [r_in * cos(θ) for θ in θ_values]
# y_inner = [r_in * sin(θ) for θ in θ_values]
# x_outer = [r_out * cos(θ) for θ in θ_values]
# y_outer = [r_out * sin(θ) for θ in θ_values]

# # Initialize the plot
# pl_disc = plot(
#     x_outer, y_outer,
#     seriestype = :shape,
#     fillcolor = :mediumorchid2,
#     linecolor = :rebeccapurple,
#     aspect_ratio = :equal,
#     xlabel = "x",
#     ylabel = "y",
#     title = "Disc Hits Visualization",
#     label = "Disc",
#     xlim = (-20, 20),
#     ylim = (-20, 20),
# )

# # Fill the inner disc area
# plot!(
#     x_inner, y_inner,
#     seriestype = :shape,
#     fillcolor = :white,
#     linecolor = :rebeccapurple,
#     label = false
# )

# # Overlay the hit points
# scatter!(
#     pl_disc,
#     x_hits, y_hits,
#     color = :tan1,
#     markerstrokecolor = :black,
#     markersize = 6,
#     label = "Disc Hits",
# )

# # Optionally, plot the geodesic paths that hit the disc
# for i in 1:length(b_values)
#     b = b_values[i]
#     # Check if this b corresponds to a disc hit
#     hit_indices = findall(x -> x[1] == b, disc_hits)
#     if !isempty(hit_indices)
#         # Plot the trajectory
#         plot!(
#             pl_disc,
#             x_vals_all[i], y_vals_all[i],
#             lw = 1.0,
#             linecolor = :lightsteelblue,
#             label = false
#         )
#     end
# end

# # Overlay the event horizon
# circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
# circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]
# plot!(
#     circle_x, circle_y,
#     lw = 2,
#     color = :purple4,
#     label = "Event horizon"
# )

# # Overlay the ISCO
# isco_x = [r_in * cos(θ) for θ in 0:0.01:2π]
# isco_y = [r_in * sin(θ) for θ in 0:0.01:2π]
# plot!(
#     isco_x, isco_y,
#     lw = 2,
#     color = :mediumslateblue,
#     label = "ISCO"
# )


# # Display the plot
# display(pl_disc)
