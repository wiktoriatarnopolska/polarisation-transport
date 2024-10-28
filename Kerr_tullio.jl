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
    g_tt = -(1 - 2 * M * r / Σ)
    g_tϕ = -2 * M * a * r * sin(θ)^2 / Σ
    g_rr = Σ / Δ
    g_θθ = Σ
    g_ϕϕ = ((r^2 + a^2)^2 - Δ * a^2 * sin(θ)^2) * sin(θ)^2 / Σ

    g = [
        g_tt    0       0       g_tϕ;
        0       g_rr    0       0;
        0       0       g_θθ    0;
        g_tϕ    0       0       g_ϕϕ
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

# Create initial velocity vector (t, r, θ, ϕ components)
v_initial = [0.0, v_r, v_θ, v_ϕ]

# Extract spatial components
v_space = v_initial[2:4]
g_space = g0[2:4, 2:4]

# Compute spatial term: v_i * g_ij * v_j
spatial_term = v_space' * g_space * v_space  # This gives a scalar

# Compute v_t using the normalization condition
v_t = sqrt((-1.0 - spatial_term) / g_tt)

# Complete the initial velocity vector
v = [v_t; v_r; v_θ; v_ϕ]

# Verify the normalization condition
norm = g0[1,1] * v[1]^2 + g0[2,2] * v[2]^2 + g0[3,3] * v[3]^2 + g0[4,4] * v[4]^2
println("Normalization check: ", norm)  # Should be close to -1

# Combine position and velocity into initial condition vector (8 elements)
# *** Correction made here: using v[1], v[2], v[3], v[4] instead of v_initial ***
u0 = [t0, r0, θ0, ϕ0, v[1], v[2], v[3], v[4]]

# Event horizon radius for the Kerr metric
r_horizon = M + sqrt(M^2 - a^2)

# Analytical Christoffel symbols for Schwarzschild metric
function compute_christoffel_analytical(r, θ)
    M = 1.0
    Γ = zeros(4, 4, 4)

    # Precompute common terms
    f = 1 - 2 * M / r
    f_inv = 1 / f

    # Non-zero Christoffel symbols
    # Γ^t_{tr} = Γ^t_{rt}
    Γ[1,2,1] = M / (r^2 * f)
    Γ[1,1,2] = Γ[1,2,1]  # Symmetry in lower indices

    # Γ^r_{tt}
    Γ[2,1,1] = f * M / r^2

    # Γ^r_{rr}
    Γ[2,2,2] = M / (r^2 * f)

    # Γ^r_{θθ}
    Γ[2,3,3] = -r * f

    # Γ^r_{ϕϕ}
    Γ[2,4,4] = -r * f * sin(θ)^2

    # Γ^θ_{rθ} = Γ^θ_{θr}
    Γ[3,2,3] = 1 / r
    Γ[3,3,2] = Γ[3,2,3]  # Symmetry

    # Γ^θ_{ϕϕ}
    Γ[3,4,4] = -sin(θ) * cos(θ)

    # Γ^ϕ_{rϕ} = Γ^ϕ_{ϕr}
    Γ[4,2,4] = 1 / r
    Γ[4,4,2] = Γ[4,2,4]  # Symmetry

    # Γ^ϕ_{θϕ} = Γ^ϕ_{ϕθ}
    Γ[4,3,4] = cos(θ) / sin(θ)
    Γ[4,4,3] = Γ[4,3,4]  # Symmetry

    return Γ
end


# Define the ODE function using an in-place form
function intprob!(du, u, p, t)
    # Current position and velocity
    x = u[1:4]
    v = u[5:8]

    r_val, θ_val = x[2], x[3]

    # Initialize du
    du .= 0.0

    # Terminate the integration if r crosses the event horizon
    if r_val <= 2 * M
        du[1:4] .= v
        return
    end

    # Compute analytical Christoffel symbols
    Γ = compute_christoffel_analytical(r_val, θ_val)

    # Compute the derivatives of velocity
    dv = zeros(4)
    for i in 1:4
        sum_ = 0.0
        for j in 1:4
            for k in 1:4
                sum_ += Γ[i,j,k] * v[j] * v[k]
            end
        end
        dv[i] = -sum_
    end

    # Assign derivatives to du
    du[1:4] .= v
    du[5:8] .= dv
end

tspan = (0.0, 1000.0)

# Solve the ODE using in-place form
prob = ODEProblem(intprob!, u0, tspan)
sol = solve(prob, Vern9(), abstol=1e-12, reltol=1e-12, dtmax=0.1)

# Extract radial and azimuthal values
r_vals = [sol[i][2] for i in 1:length(sol)]
ϕ_vals = [sol[i][4] for i in 1:length(sol)]

# Plot the geodesic path
pl = plot(
    θ -> 2 * M,  # Event horizon radius for Schwarzschild
    0:0.01:2π,
    lw = 2,
    legend = true,
    proj = :polar,
    color = :black,
    ylim=(0.0, 15.0),
    label = "Event horizon"
)

plot!(ϕ_vals, r_vals, lw=2, label="Geodesic Path", color=:blue)
display(pl)
