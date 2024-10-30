using LinearAlgebra, DifferentialEquations, Plots
using Roots

# Constants
M = 1.0
a = 0.0  # Set this to a non-zero value for Kerr spacetime (rotating black hole)

# Initial conditions
r0 = 3.0
θ0 = π / 2
ϕ0 = 0.0
λ0 = 0.0

# Kerr Metric components as a function of r and θ
function metric(r, θ)
    Δ = r^2 - 2 * r + a^2
    Σ = r^2 + a^2 * cos(θ)^2

    g_tt = -(1 - 2 * M * r / Σ)
    g_rr = Σ / Δ
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


# Set initial velocity components (adjust these values as desired)
v_r = 0.0    # Radial velocity (negative for inward motion)
v_θ = 0.0     # Polar velocity
v_ϕ = 0.3     # Azimuthal (angular) velocity

# Extract metric components
g_tt = g0[1,1]
g_tϕ = g0[1,4]
g_rr = g0[2,2]
g_θθ = g0[3,3]
g_ϕϕ = g0[4,4]

# Compute coefficients for quadratic equation in v_t
A = g_tt
B = g_tϕ * v_ϕ
C = g_ϕϕ * v_ϕ^2 + g_rr * v_r^2 + g_θθ * v_θ^2

# Solve quadratic equation A * v_t^2 + 2 * B * v_t + C = 0
Δ = (2 * B)^2 - 4 * A * C
if Δ < 0
    error("No real solution for v_t. Check initial conditions.")
end

v_t = (-2 * B - sqrt(Δ)) / (2 * A)  # Choose the negative root for future-directed motion

# Complete the initial velocity vector
v = [v_t; v_r; v_θ; v_ϕ]

# Verify the normalization condition
norm_ = g_tt * v[1]^2 + g_rr * v[2]^2 + g_θθ * v[3]^2 + g_ϕϕ * v[4]^2 +
g_tϕ * v[1] * v[4]
println("Null normalization check: ", norm_)  # Should be close to 0

# Combine position and velocity into initial condition vector (8 elements)
u0 = [λ0, r0, θ0, ϕ0, v[1], v[2], v[3], v[4]]

# Event horizon for the Kerr metric
r_horizon = sqrt(1 - a^2)

# Compute the Kerr Christoffel symbols analytically
function compute_christoffel_analytical(r, θ)
    Δ = r^2 - 2 * r + a^2
    Σ = r^2 + a^2 * cos(θ)^2
    Γ = zeros(4, 4, 4)
    f = r^2 - a^2 * cos(θ)^2
    A = (r^2 + a^2) * Σ + 2 * a^2 * r * sin(θ)^2

    Γ[1, 1, 2] = (r^2 + a^2) * f / Σ^2 * Δ
    Γ[1, 2, 1] = Γ[1, 1, 2]
    Γ[1, 1, 3] = - 2 * a^2* r * sin(θ) * cos(θ) / Σ^2
    Γ[1, 3, 1] = Γ[1, 1, 3]

    Γ[2, 1, 1] = Δ * f / Σ^3
    Γ[2, 1, 4] = - Δ * a * sin(θ)^2 * f / Σ^3
    Γ[2, 4, 1] = Γ[2, 1, 4]
    Γ[2, 2, 2] = (r * a^2 * sin(θ)^2 - f) / Σ * Δ
    Γ[2, 2, 3] = - a^2 * sin(θ) * cos(θ) / Σ
    Γ[2, 3, 2] = Γ[2, 2, 3]
    Γ[2, 3, 3] = - r * Δ / Σ

    Γ[4, 3, 4] = ((cos(θ) / sin(θ)) / Σ^2) * (Σ^2 + 2 * a^2 * r * sin(θ)^2)
    Γ[4, 4, 3] = Γ[4, 3, 4]
    
    Γ[3, 1, 1] = - 2 * a^2 * r * sin(θ) * cos(θ) / Σ^3

    Γ[4, 1, 2] = a * f / (Σ^2 * Δ)
    Γ[4, 2, 1] = Γ[4, 1, 2]

    Γ[4, 1, 3] = - 2 * a * r * (cos(θ) / sin(θ)) / Σ^2
    Γ[4, 3, 1] = Γ[4, 1, 3]

    Γ[3, 1, 4] = 2 * a * r * (r^2 + a^2) * sin(θ) * cos(θ) / Σ^3
    Γ[3, 4, 1] = Γ[3, 1, 4]

    Γ[3, 2, 2] = a^2 * sin(θ) * cos(θ) / Σ * Δ

    Γ[3, 2, 3] = r / Σ
    Γ[3, 3, 2] = Γ[3, 2, 3]

    Γ[3, 3, 3] = - a^2 * sin(θ) * cos(θ) / Σ

    Γ[1, 3, 4] = 2 * a^3 * r * sin(θ)^3 * cos(θ) / Σ^2
    Γ[1, 4, 3] = Γ[1, 3, 4]

    Γ[1, 2, 4] = a * sin(θ)^2 * (a^2 * cos(θ)^2 *(a^2 - r^2) - r^2 * (a^2 + 3 * r^2)) / (2 * Σ^2 * Δ)
    Γ[1, 4, 2] = Γ[1, 2, 4]

    Γ[4, 2, 4] = (r * Σ^2 + (a^4 * sin(θ)^2 * cos(θ)^2 - r^2(Σ + r^2 + a^2))) / (Σ * Δ)
    Γ[4, 4, 2] = Γ[4, 2, 4]

    Γ[2, 4, 4] = (Δ * sin(θ)^2 / (Σ^3)) * (- r * Σ^2 + a^2 * sin(θ)^2 * f)
    
    Γ[3, 4, 4] = (- sin(θ) * cos(θ) / Σ^3) * (A * Σ + (r^2 + a^2) * 2 * a^2 * r * sin(θ)^2 )

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

# Set up and solve the ODE problem
tspan = (0.0, 1000.0)
prob = ODEProblem(intprob!, u0, tspan)
sol = solve(prob, Tsit5(), abstol=1e-14, reltol=1e-14, dtmax=0.01)

# # Initialize arrays to store conserved quantities
# E_vals = []
# L_vals = []

for i in 1:length(sol)
    x = sol[i][1:4]
    v = sol[i][5:8]

    r = x[2]
    θ = x[3]

    # Compute the metric components at the current position
    g = metric(r, θ)
    g_tt = g[1,1]
    g_ϕϕ = g[4,4]

    # # Compute conserved quantities
    # E = -g_tt * v[1]  # Energy-like quantity
    # L = g_ϕϕ * v[4]   # Angular momentum-like quantity

    # push!(E_vals, E)
    # push!(L_vals, L)
end

# Extract radial and azimuthal values for plotting
r_vals = [sol[i][2] for i in 1:length(sol)]
ϕ_vals = [sol[i][4] for i in 1:length(sol)]

# Plot the geodesic path and event horizon for Kerr metric
pl = plot(
    θ -> r_horizon,
    0:0.01:2π,
    lw = 2,
    legend = true,
    proj = :polar,
    color = :black,
    ylim=(0.0, 10),
    label="Event Horizon"
)

plot!(ϕ_vals, r_vals, lw=2, label="Photon Path", color=:blue)
