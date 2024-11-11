using LinearAlgebra, DifferentialEquations, Plots

# Constants
M = 1.0  # Mass of the black hole
a = 0.9

observer = (1000.0, deg2rad(60), 0.0, 0.0)  # θ_obs = 60 degrees

# Event horizon radius
r_horizon = 1 + sqrt(1 - a^2)

# Kerr Metric components as a function of r and θ
function metric(r, θ)
    Δ = r^2 - 2 * r + a^2
    Σ = r^2 + a^2 * cos(θ)^2
    f = 1 - 2 * M / r

    g_tt = -(f * r^2 / Σ)
    g_rr = 1/f
    g_θθ = Σ
    g_ϕϕ = (r^2 + a^2 + ((2 * M * r * a^2 * sin(θ)^2 )/ Σ)) * sin(θ)^2
    g_tϕ = -2 * M * r * a * sin(θ)^2 / Σ

    g = [
        g_tt     0        0        g_tϕ;
        0        g_rr     0        0;
        0        0        g_θθ     0;
        g_tϕ     0        0        g_ϕϕ
    ]
    return g
end

# Compute the Kerr Christoffel symbols analytically
function compute_christoffel_analytical(r, θ)
    Δ = r^2 - 2 * r + a^2
    Σ = r^2 + a^2 * cos(θ)^2
    Γ = zeros(4, 4, 4)
    #f = r^2 - a^2 * cos(θ)^2
    #f = 1 - 2 / r

    A = (r^2 + a^2) * Σ + 2 * a^2 * r * sin(θ)^2

    Γ[1, 1, 2] = (r^2 + a^2) * (r^2 - a^2 * cos(θ)^2) / ((Σ^2) * Δ)
    Γ[1, 2, 1] = Γ[1, 1, 2]

    Γ[1, 1, 3] = - 2 * a^2* r * sin(θ) * cos(θ) / Σ^2
    Γ[1, 3, 1] = Γ[1, 1, 3]

    Γ[2,1,1] = Δ * (r^2 - a^2 * cos(θ)^2) / (Σ^3)

    Γ[2, 1, 4] = - Δ * a * sin(θ)^2 * (r^2 - a^2 * cos(θ)^2) / Σ^3
    Γ[2, 4, 1] = Γ[2, 1, 4]

    Γ[2, 2, 2] = (r * a^2 * sin(θ)^2 - (r^2 - a^2 * cos(θ)^2)) / (Σ * Δ)

    Γ[2, 2, 3] = - a^2 * sin(θ) * cos(θ) / Σ
    Γ[2, 3, 2] = Γ[2, 2, 3]

    Γ[2, 3, 3] = - r * Δ / Σ

    Γ[2, 4, 4] = (Δ * sin(θ)^2 / (Σ^3)) * (- r * Σ^2 + a^2 * sin(θ)^2 * (r^2 - a^2 * cos(θ)^2))

    Γ[3, 2, 3] = r / Σ
    Γ[3, 3, 2] = Γ[3, 2, 3]

    Γ[3, 4, 4] = (- sin(θ) * cos(θ) / (Σ^3)) * (A * Σ + (r^2 + a^2) * 2 * a^2 * r * sin(θ)^2 )

    Γ[4, 2, 4] = (r * Σ^2 + (a^4 * sin(θ)^2 * cos(θ)^2 - r^2*(Σ + r^2 + a^2))) / (Σ^2 * Δ)
    Γ[4, 4, 2] = Γ[4, 2, 4]

    Γ[4, 3, 4] = ((cos(θ) / sin(θ)) / Σ^2) * (Σ^2 + 2 * a^2 * r * sin(θ)^2)
    Γ[4, 4, 3] = Γ[4, 3, 4]

    Γ[3, 1, 1] = - 2 * a^2 * r * sin(θ) * cos(θ) / Σ^3

    Γ[4, 1, 2] = a * (r^2 - a^2 * cos(θ)^2) / (Σ^2 * Δ)
    Γ[4, 2, 1] = Γ[4, 1, 2]

    Γ[4, 1, 3] = - 2 * a * r * (cos(θ) / sin(θ)) / Σ^2
    Γ[4, 3, 1] = Γ[4, 1, 3]

    Γ[3, 1, 4] = 2 * a * r * (r^2 + a^2) * sin(θ) * cos(θ) / Σ^3
    Γ[3, 4, 1] = Γ[3, 1, 4]

    Γ[3, 2, 2] = a^2 * sin(θ) * cos(θ) / (Σ * Δ)

    Γ[3, 3, 3] = - a^2 * sin(θ) * cos(θ) / Σ

    Γ[1, 3, 4] = 2 * a^3 * r * sin(θ)^3 * cos(θ) / Σ^2
    Γ[1, 4, 3] = Γ[1, 3, 4]

    Γ[1, 2, 4] = a * sin(θ)^2 * (a^2 * cos(θ)^2 * (a^2 - r^2) - r^2 * (a^2 + 3 * r^2)) / (Σ^2 * Δ)
    Γ[1, 4, 2] = Γ[1, 2, 4]

    return Γ
end

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

tspan = (0.0, 5000.0)

# Define a range of impact parameters
b_values = collect(-10.0:1:10.0)

# Arrays to store trajectories
x_vals_all = []
y_vals_all = []

# Disc parameters
r_in = isco_radius(0.9)       # Inner radius of the disc
r_out = 10.0                  # Outer radius of the disc
N_r = 50                      # Radial grid points
N_φ = 50                      # Azimuthal grid points

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

    # Solve the ODE
    prob = ODEProblem(intprob!, u0, tspan)
    sol = solve(prob, Tsit5(), abstol=1e-12, reltol=1e-12, dtmax=1.0)

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
for i in 1:min(5, length(disc_hits))
    println(disc_hits[i])
end

# Now you can use disc_hits to analyze or visualize the geodesics that hit the disc
for (b, r_hit, ϕ_hit) in disc_hits
    println("Geodesic with impact parameter b = $b hits the disc at (r, ϕ) = ($r_hit, $ϕ_hit)")
end


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
    fillcolor = :lightblue,
    linecolor = :blue,
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    title = "Disc Hits Visualization",
    legend = false,
    xlim = (-20, 20),
    ylim = (-20, 20),
)

# Fill the inner disc area
plot!(
    x_inner, y_inner,
    seriestype = :shape,
    fillcolor = :white,
    linecolor = :blue
)

# Overlay the hit points
scatter!(
    pl_disc,
    x_hits, y_hits,
    color = :red,
    markerstrokecolor = :black,
    markersize = 6,
    label = "Disc Hits",
)

# Optionally, plot the geodesic paths that hit the disc
for i in 1:length(b_values)
    b = b_values[i]
    # Check if this b corresponds to a disc hit
    hit_indices = findall(x -> x[1] == b, disc_hits)
    if !isempty(hit_indices)
        # Plot the trajectory
        plot!(
            pl_disc,
            x_vals_all[i], y_vals_all[i],
            lw = 1.0,
            linecolor = :gray,
        )
    end
end

# Display the plot
display(pl_disc)
 