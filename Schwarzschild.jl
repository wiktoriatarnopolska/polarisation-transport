using LinearAlgebra, DifferentialEquations, Plots

# Constants
M = 1.0
a = 0.0  # Schwarzschild spacetime (non-rotating black hole)

# Initial conditions
r0 = 10.0       # Initial radial distance
θ0 = π / 2      # Equatorial plane
ϕ0 = 0.0        # Initial azimuthal angle
t0 = 0.0        # Initial time

# Calculate specific energy E and angular momentum L for circular orbit
E = (1 - 2 * M / r0) / sqrt(1 - 3 * M / r0)
L = sqrt(M * r0) / sqrt(1 - 3 * M / r0)

# Metric components as a function of r and θ
function metric(r, θ)
    Σ = r^2 + a^2 * cos(θ)^2
    Δ = r^2 - 2 * M * r + a^2

    # Metric components
    g_tt = -(1 - 2 * M / r)
    g_rr = 1 / (1 - 2 * M / r)
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

# Compute the metric at the initial position
g0 = metric(r0, θ0)

# Compute the initial 4-velocity components
g_tt = g0[1,1]
g_ϕϕ = g0[4,4]

v_r = 0.0
v_θ = 0.0
v_ϕ = L / g_ϕϕ  # v^ϕ = L / g_ϕϕ

# Compute spatial term: v_i * g_ij * v_j
spatial_term = g_ϕϕ * v_ϕ^2  # Since v_r = v_θ = 0

# Compute v_t using the normalization condition
v_t = sqrt((-1.0 - spatial_term) / g_tt)

# Complete the initial velocity vector
v = [v_t; v_r; v_θ; v_ϕ]

# Verify the normalization condition
norm = g0[1,1] * v[1]^2 + g0[4,4] * v[4]^2
println("Normalization check: ", norm)  # Should be close to -1

# Combine position and velocity into initial condition vector (8 elements)
u0 = [t0, r0, θ0, ϕ0, v[1], v[2], v[3], v[4]]

# Event horizon radius for the Schwarzschild metric
r_horizon = 2 * M

# Analytical Christoffel symbols for Schwarzschild metric
function compute_christoffel_analytical(r, θ)
    M = 1.0
    Γ = zeros(4, 4, 4)

    # Precompute common terms
    f = 1 - 2 * M / r

    # Non-zero Christoffel symbols
    Γ[1,2,1] = M / (r^2 * f)
    Γ[1,1,2] = Γ[1,2,1]  # Symmetry in lower indices

    Γ[2,1,1] = f * M / r^2

    Γ[2,2,2] = M / (r^2 * f)

    Γ[2,3,3] = -r * f

    Γ[2,4,4] = -r * f * sin(θ)^2

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
function intprob!(du, u, p, τ)
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

tspan = (0.0, 1000.0)

# Solve the ODE using in-place form
prob = ODEProblem(intprob!, u0, tspan)
sol = solve(prob, Vern9(), abstol=1e-14, reltol=1e-14, dtmax=0.01)

# Extract radial and azimuthal values
r_vals = [sol[i][2] for i in 1:length(sol)]
ϕ_vals = [sol[i][4] for i in 1:length(sol)]

# Plot the geodesic path
pl = plot(
    θ -> r_horizon,  # Event horizon radius for Schwarzschild
    0:0.01:2π,
    lw = 2,
    legend = true,
    proj = :polar,
    color = :black,
    ylim=(0.0, 15.0),
    label="Event Horizon"
)

plot!(ϕ_vals, r_vals, lw=2, label="Geodesic Path", color=:blue)
display(pl)

