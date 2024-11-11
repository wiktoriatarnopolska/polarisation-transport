using LinearAlgebra, DifferentialEquations, Plots

# Constants
M = 1.0  # Mass of the black hole
a = 0.9

observer = (1000.0, deg2rad(90), 0.0, 0.0)

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

λ_vals_all = []
E_vals_all = []
L_vals_all = []
Q_vals_all = []

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
                               # this ensures that the photon starts with 
                               # an angular momentum corresponding to the 
                               # impact parameter b at the initial large radius.

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

    # Extract positions
    r_vals = [sol[i][2] for i in 1:length(sol)]
    θ_vals = [sol[i][3] for i in 1:length(sol)]
    ϕ_vals = [sol[i][4] for i in 1:length(sol)]

    # Convert to Cartesian coordinates
    x_vals = [r * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
    y_vals = [r * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]

    # Append trajectories
    push!(x_vals_all, x_vals)
    push!(y_vals_all, y_vals)

    # Extract conserved quantities for each step
    λ_vals = sol.t
    E_vals = Float64[]
    L_vals = Float64[]
    Q_vals = Float64[]
    
    for i in 1:length(sol)
        x = sol[i][1:4]  # Position
        v = sol[i][5:8]  # Velocity
        r, θ, ϕ = x[2], x[3], x[4]

        g = metric(r, θ)  # Metric at current position
        g_tt, g_tϕ, g_rr, g_θθ, g_ϕϕ = g[1, 1], g[1, 4], g[2, 2], g[3, 3], g[4, 4]

        # Compute conserved quantities
        E = -g_tt * v[1] - g_tϕ * v[4]
        L = g_tϕ * v[1] + g_ϕϕ * v[4]
        Q = (g_θθ * v[3])^2 + cos(θ)^2 * ( - a^2 * (g_tt * v[1])^2 + (g_ϕϕ * v[4] + g_tϕ * v[1])^2 / sin(θ)^2)

        push!(E_vals, E)
        push!(L_vals, L)
        push!(Q_vals, Q)
    end

    # Store the affine parameter and conserved quantities for this trajectory
    push!(λ_vals_all, λ_vals)
    push!(E_vals_all, E_vals)
    push!(L_vals_all, L_vals)
    push!(Q_vals_all, Q_vals)

end

# Initialize the plot
pl_cartesian = plot(
    xlabel = "x",
    ylabel = "y",
    aspect_ratio = :equal,
    title = "Kerr metric a = $a",
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

pl_conserved = plot(
    xlabel="Affine Parameter λ",
    ylabel="Quantity",
    title="Q, E and L conservation",
    legend=:outerright,
    legendfontsize=3
)

# Loop through each trajectory and plot
for i in 1:length(b_values)
    λ_vals = λ_vals_all[i]
    E_vals = E_vals_all[i]
    L_vals = L_vals_all[i]
    Q_vals = Q_vals_all[i]

    # Plot Energy-like quantity E(λ)
    plot!(
        pl_conserved,
        λ_vals, E_vals,
        label="E(λ) for b=$(b_values[i])",
        lw=1.5,
        color=:blue
    )

    # Plot Angular Momentum-like quantity L(λ)
    plot!(
        pl_conserved,
        λ_vals, L_vals,
        label="L(λ) for b=$(b_values[i])",
        lw=1.5,
        color=:red
    )

    # Plot Carter Constant Q(λ)
    plot!(
        pl_conserved,
        λ_vals, Q_vals,
        label="Q(λ) for b=$(b_values[i])",
        lw=1.5,
        color=:green
    )
end

# Display the plot
display(pl_conserved)
