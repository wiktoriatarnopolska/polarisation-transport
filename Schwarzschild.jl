using LinearAlgebra, DifferentialEquations, Plots

# Constants
M = 1.0
a = 0.0  # Schwarzschild spacetime (non-rotating black hole)

# Initial conditions
#r0 = 10.0       # Initial radial distance
r0 = 1000.0
θ0 = π / 2      # Equatorial plane
ϕ0 = 0.0        # Initial azimuthal angle
λ0 = 0.0        # Initial affine parameter

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

# Set initial velocity components (adjust these values as desired)
v_r = 0.0    # Radial velocity (negative for inward motion)
v_θ = 0.0     # Polar velocity
v_ϕ = 0.3     # Azimuthal (angular) velocity

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

# Verify the normalization condition
norm = g_tt * v[1]^2 + g_rr * v[2]^2 + g_θθ * v[3]^2 + g_ϕϕ * v[4]^2
println("Null normalization check: ", norm)  # Should be close to 0

# Combine position and velocity into initial condition vector (8 elements)
u0 = [λ0, r0, θ0, ϕ0, v[1], v[2], v[3], v[4]]

# Event horizon radius for the Schwarzschild metric
r_horizon = 2

# Analytical Christoffel symbols for Schwarzschild metric
function compute_christoffel_analytical(r, θ)
    M = 1.0
    Γ = zeros(4, 4, 4)

    # Precompute common terms
    f = 1 - 2 * M / r

    # Non-zero Christoffel symbols
    Γ[1,2,1] = 1 / (r^2 * f)
    Γ[1,1,2] = Γ[1,2,1]  # Symmetry in lower indices

    Γ[2,1,1] = f / r^2

    Γ[2,2,2] = 1 / (r^2 * f)

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

tspan = (0.0, 1000.0)

# Solve the ODE using in-place form
prob = ODEProblem(intprob!, u0, tspan)
sol = solve(prob, Tsit5(), abstol=1e-14, reltol=1e-14, dtmax=0.01)

# Initialize arrays to store conserved quantities
E_vals = []
L_vals = []

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
    ylim=(0.0, r0 + 1.0),
    label="Event Horizon"
)

plot!(ϕ_vals, r_vals, lw=2, label="Photon Path", color=:blue)


# # Extract affine parameter values
# λ_vals = sol.t

# # Plot Energy-like and Momentum-like Quantity
# plot(λ_vals, L_val1s, xlabel="Affine Parameter λ", ylabel="Angular Momentum-like Quantity L", label="L(λ)", colour =:blue )
# plot(λ_vals, E_vals, xlabel="Affine Parameter λ", ylabel="Energy-like Quantity E", label="E(λ)", legend=:bottomright, colour =:red)
# plot!(λ_vals, L_vals, title = "Conservation of E and L", xlabel="Affine Parameter λ", ylabel="Angular Momentum-like Quantity L", label="L(λ)", colour =:blue )

# x_vals = [r * cos(ϕ) * sin(θ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
# y_vals = [r * sin(ϕ) * sin(θ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
# z_vals = [r * cos(θ) for (r, θ) in zip(r_vals, θ_vals)]

# # Plot the photon path
# plot(x_vals, y_vals, z_vals, lw=2, label="Photon Path", aspect_ratio=:equal)
# plot!(xlabel="x", ylabel="y", zlabel="z", title="Photon Orbit in Kerr Spacetime")
