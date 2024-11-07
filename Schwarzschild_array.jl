using LinearAlgebra, DifferentialEquations, Plots

# Constants
M = 1.0  # Mass of the black hole

# Event horizon radius
r_horizon = 2 * M  # Schwarzschild radius

# Metric components as a function of r and θ
function metric(r, θ)
    f = 1 - 2 * M / r
    g_tt = -f
    g_rr = 1 / f
    g_θθ = r^2
    g_ϕϕ = r^2 * sin(θ)^2

    g = [
        g_tt    0       0       0;
        0       g_rr    0       0;
        0       0       g_θθ    0;
        0       0       0       g_ϕϕ
    ]
    return g
end

# Analytical Christoffel symbols for Schwarzschild metric
function compute_christoffel_analytical(r, θ)

    Γ = zeros(4, 4, 4)

    # Non-zero Christoffel symbols
    Γ[1,2,1] = 1 / (r * (r - 2))
    Γ[1,1,2] = Γ[1,2,1]  # Symmetry in lower indices

    Γ[2,1,1] = (1 - 2 / r) / r^2

    Γ[2,2,2] = - 1 / (r^2 * (1 - 2 / r))

    Γ[2,3,3] = -r * (1 - 2 / r)

    Γ[2,4,4] = -r * (1 - 2 / r) * sin(θ)^2

    Γ[3,2,3] = 1 / r
    Γ[3,3,2] = Γ[3,2,3]  # Symmetry

    Γ[3,4,4] = -sin(θ) * cos(θ)

    Γ[4,2,4] = 1 / r
    Γ[4,4,2] = Γ[4,2,4]  # Symmetry

    Γ[4,3,4] = cos(θ) / sin(θ)
    Γ[4,4,3] = Γ[4,3,4]  # Symmetry

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

# Observer position at large negative x
x_observer = -1000.0  # Negative x-axis
y_observer = 0.0
z_observer = 0.0  # Equatorial plane

# Define a range of impact parameters
b_values = collect(-10.0:1:10.0)

# Arrays to store trajectories
x_vals_all = []
y_vals_all = []

for b in b_values
    # Compute initial position in Cartesian coordinates
    x0_cartesian = [x_observer, y_observer, z_observer]

    # Compute initial direction in Cartesian coordinates
    α = atan(b / abs(x_observer))  # Angle with respect to x-axis
    v_x = cos(α)
    v_y = sin(α)
    v_z = 0.0  # Equatorial plane

    # Since x_observer is negative, photons move towards increasing x
    v_x = abs(v_x)

    # Normalize the spatial velocity components
    v_magnitude = sqrt(v_x^2 + v_y^2 + v_z^2)
    v_x /= v_magnitude
    v_y /= v_magnitude

    # Convert initial position to spherical coordinates
    r0 = sqrt(x_observer^2 + y_observer^2 + z_observer^2)
    θ0 = acos(z_observer / r0)
    ϕ0 = atan(y_observer, x_observer)

    # Convert initial velocities to spherical coordinates
    # Compute basis vectors at initial position
    sin_θ0 = sin(θ0)
    cos_θ0 = cos(θ0)
    sin_ϕ0 = sin(ϕ0)
    cos_ϕ0 = cos(ϕ0)

    e_r = [sin_θ0 * cos_ϕ0, sin_θ0 * sin_ϕ0, cos_θ0]
    e_θ = [cos_θ0 * cos_ϕ0, cos_θ0 * sin_ϕ0, -sin_θ0]
    e_ϕ = [-sin_ϕ0, cos_ϕ0, 0.0]

    # Compute v_r, v_θ, v_ϕ
    v_cartesian = [v_x, v_y, v_z]

    v_r = dot(v_cartesian, e_r)
    v_θ = dot(v_cartesian, e_θ) / r0
    v_ϕ = dot(v_cartesian, e_ϕ) / (r0 * sin_θ0)

    # Compute the metric at the initial position
    g0 = metric(r0, θ0)
    g_tt = g0[1,1]
    g_rr = g0[2,2]
    g_θθ = g0[3,3]
    g_ϕϕ = g0[4,4]

    # Compute spatial term: v_i * g_ij * v_j
    spatial_term = g_rr * v_r^2 + g_θθ * v_θ^2 + g_ϕϕ * v_ϕ^2

    # Compute v_t using the null normalization condition
    v_t = sqrt(spatial_term / (-g_tt))

    # Complete the initial velocity vector
    v = [v_t; v_r; v_θ; v_ϕ]

    # Initial affine parameter
    λ0 = 0.0

    # Combine position and velocity into initial condition vector
    u0 = [λ0, r0, θ0, ϕ0, v[1], v[2], v[3], v[4]]

    # Solve the ODE
    prob = ODEProblem(intprob!, u0, tspan)
    sol = solve(prob, Tsit5(), abstol=1e-12, reltol=1e-12, dtmax=1.0)

    # Extract positions
    r_vals = [sol[i][2] for i in 1:length(sol)]
    θ_vals = [sol[i][3] for i in 1:length(sol)]
    ϕ_vals = [sol[i][4] for i in 1:length(sol)]

    # Convert to Cartesian coordinates
    x_vals = [r * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
    y_vals = [r * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]

    # Append trajectories to the arrays
    push!(x_vals_all, x_vals)
    push!(y_vals_all, y_vals)
end

# Initialize the plot
pl_cartesian = plot(
    xlabel = "x",
    ylabel = "y",
    aspect_ratio = :equal,
    title = "Schwarzschild metric a = 0",
    xlim = (-20, 20),
    ylim = (-20, 20),
    legend = false
)

# Plot each photon's path
for i in 1:length(b_values)
    plot!(
        pl_cartesian,
        x_vals_all[i], y_vals_all[i],
        lw = 1.5,
    )
end

# Overlay the event horizon
circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]
plot!(
    pl_cartesian,
    circle_x, circle_y,
    lw = 2,
    color = :black
)

# Display the plot
display(pl_cartesian)
