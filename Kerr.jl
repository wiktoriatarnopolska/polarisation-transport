using LinearAlgebra, DifferentialEquations, Plots


# Constants
M = 1.0
a = 0.9  # Set this to a non-zero value for Kerr spacetime (rotating black hole)

observer = (1000.0, deg2rad(90), 0.0, 0.0)

r0 = observer[1]
θ0 = observer[2]
ϕ0 = observer[3]
λ0 = observer[4]

# Example usage for transforming and solving for initial conditions
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
tspan = (0.0, 5000.0)
prob = ODEProblem(intprob!, u0, tspan)
sol = @time solve(prob, Tsit5(), abstol=1e-14, reltol=1e-14, dtmax=0.01)

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
    g_tϕ = g[1,4]
    g_rr = g[2,2]
    g_θθ = g[3,3]
    g_ϕϕ = g[4,4]


    # Compute conserved quantities
    E = -g_tt * v[1] - g_tϕ * v[4] # Energy-like quantity
    L = g_tϕ * v[1] + g_ϕϕ * v[4]  # Angular momentum-like quantity
    Q = (g_θθ * v[3])^2 + cos(θ)^2 * ( - a^2 * (g_tt * v[1])^2 + (g_ϕϕ * v[4] + g_tϕ * v[1])^2 / sin(θ)^2)

    push!(E_vals, E)
    push!(L_vals, L)
    push!(Q_vals, Q)
end


# Extract radial and azimuthal values for plotting
r_vals = [sol[i][2] for i in 1:length(sol)]
ϕ_vals = [sol[i][4] for i in 1:length(sol)]

# Convert polar coordinates to Cartesian coordinates
x_vals = [r * cos(ϕ) for (r, ϕ) in zip(r_vals, ϕ_vals)]
y_vals = [r * sin(ϕ) for (r, ϕ) in zip(r_vals, ϕ_vals)]

# Plot the photon path in Cartesian coordinates
pl_cartesian = plot(
    x_vals, y_vals,
    lw = 2,
    label = "Photon Path",
    color = :blue,
    xlabel = "x",
    ylabel = "y",
    aspect_ratio = :equal,
    title = "Photon Path in Cartesian Coordinates",
    xlim = (-10, 10),
    ylim = (-10, 10)
)

# Overlay event horizon in Cartesian coordinates
r_horizon = 1 + sqrt(1 - a^2)  # Event horizon radius for Kerr black hole
circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]

plot!(
    pl_cartesian,
    circle_x, circle_y,
    lw = 2,
    label = "Event Horizon",
    color = :black
)

# Display the Cartesian plot
display(pl_cartesian)


# Extract affine parameter values
λ_vals = sol.t

pl_Econserved = plot(
    xlabel="Affine Parameter λ",
    ylabel="|ΔE| (log scale)",
    title="Energy Difference |E₀ - E(λ)|",
    legend=:outerright,
    yscale=:log10,
)

abs_E_diff = abs.(E_vals .- E_vals[1])
abs_E_diff[abs_E_diff .== 0] .= 1e-10  # Replace zeros with small values

plot!(
    pl_Econserved,
    λ_vals, abs_E_diff,
    label="|ΔE|",
    lw=1.5
)


pl_Lconserved = plot(
    xlabel="Affine Parameter λ",
    ylabel="|ΔL| (log scale)",
    title="Momentum Difference |L₀ - L(λ)|",
    legend=:outerright,
    yscale=:log10,
)

abs_L_diff = abs.(L_vals .- L_vals[1])
abs_L_diff[abs_L_diff .== 0] .= 1e-10  # Replace zeros with small values

plot!(
    pl_Lconserved,
    λ_vals, abs_L_diff,
    label="|ΔL|",
    lw=1.5
)


pl_Qconserved = plot(
    xlabel="Affine Parameter λ",
    ylabel="|ΔQ| (log scale)",
    title="Carter Const. Difference |Q₀ - Q(λ)|",
    legend=:outerright,
    yscale=:log10,
)


abs_Q_diff = abs.(Q_vals .- Q_vals[1])
abs_Q_diff[abs_Q_diff .== 0] .= 1e-10  # Replace zeros with small values

plot!(
    pl_Qconserved,
    λ_vals, abs_Q_diff,
    label="|ΔQ|",
    lw=1.5
)

