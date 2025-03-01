using LinearAlgebra, DifferentialEquations, Plots
using LaTeXStrings


# Constants
M = 1.0
a = 0.0  # Schwarzschild spacetime (non-rotating black hole)

# Initial conditions
r0 = 3.0       # Initial radial distance
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

# Event horizon radius 
r_horizon = 1 + sqrt(1 - a^2)

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

tspan = (0.0, 1000.0)

# Solve the ODE using in-place form
prob = ODEProblem(intprob!, u0, tspan)
sol = solve(prob, Tsit5(), abstol=1e-12, reltol=1e-12, dtmax=0.01)

# Initialize arrays to store conserved quantities
E_vals = []
L_vals = []
Q_vals = []

for i in 1:length(sol)
    x = sol[i][1:4]
    v = sol[i][5:8]

    r = x[2]
    θ = x[3]

    # Compute the metric components at the current position
    g = metric(r, θ)
    g_tt = g[1,1]
    g_rr = g[2,2]
    g_θθ = g[3,3]
    g_ϕϕ = g[4,4]



    # Compute conserved quantities
    E = -g_tt * v[1]  # Energy-like quantity
    L = g_ϕϕ * v[4]   # Angular momentum-like quantity
    Q = (g_θθ * v[3])^2 + cos(θ)^2 * (((g_ϕϕ * v[4])^2 / sin(θ)^2))


    push!(E_vals, E)
    push!(L_vals, L)
    push!(Q_vals, Q)
end


# Extract radial and azimuthal values
r_vals = [sol[i][2] for i in 1:length(sol)]
ϕ_vals = [sol[i][4] for i in 1:length(sol)]

# First plot the photon path, then overlay the event horizon
pl = plot(
    ϕ_vals, r_vals,
    lw = 2,
    label = "Photon Path",
    color = :blue,
    proj = :polar,
    ylim = (0.0, 15),
)

plot!(
    θ -> r_horizon,
    0:0.01:2π,
    lw = 2,
    label = "Event Horizon",
    color = :black
)

# Convert polar coordinates to Cartesian coordinates
x_vals = [r * cos(ϕ) for (r, ϕ) in zip(r_vals, ϕ_vals)]
y_vals = [r * sin(ϕ) for (r, ϕ) in zip(r_vals, ϕ_vals)]

# Plot the photon path in Cartesian coordinates
pl_cartesian = plot(
    x_vals, y_vals,
    lw = 2,
    label = L"\textrm{Photon \, \, path}",
    color = :indigo,
    xlabel = L"x",
    ylabel = L"y",
    aspect_ratio=:equal,
    xlim = (-5, 5),
    ylim = (-5, 5)
)

# Overlay event horizon in Cartesian coordinates
circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]

plot!(
    pl_cartesian,
    circle_x, circle_y,
    lw = 2,
    label = L"\textrm{Event \, horizon}",
    color = :black
)

# Display both plots
plot(pl, pl_cartesian, layout = (1, 2), size = (1200, 600))

λ_vals = sol.t

abs_E_diff = abs.(E_vals .- E_vals[1])
abs_E_diff[abs_E_diff .== 0] .= 1e-10  # Replace zeros with small values

abs_L_diff = abs.(L_vals .- L_vals[1])
abs_L_diff[abs_L_diff .== 0] .= 1e-10  # Replace zeros with small values

abs_Q_diff = abs.(Q_vals .- Q_vals[1])
abs_Q_diff[abs_Q_diff .== 0] .= 1e-10  # Replace zeros with small values

pl_conserved = plot(
    xlabel=L"\textrm{Affine \, parameter \,} λ",
    ylabel=L"|Δ\textrm{Quantity}| \, \, (\textrm{log \, scale})",
    legend=:outerright,
    yscale=:log10,
)

plot!(
    pl_conserved,
    λ_vals, abs_E_diff,
    label=L"|Δ\textrm{E}|",
    lw=2,
    color=:firebrick
)


plot!(
    pl_conserved,
    λ_vals, abs_L_diff,
    label=L"|Δ\textrm{L}|",
    lw=2,
    linestyle=:dash,
    color=:blue2
)

plot!(
    pl_conserved,
    λ_vals, abs_Q_diff,
    label=L"|Δ\textrm{Q}|",
    lw=2,
    linestyle=:dot,
    color=:seagreen
)

plot(pl_cartesian, pl_conserved, layout = (1, 2))
