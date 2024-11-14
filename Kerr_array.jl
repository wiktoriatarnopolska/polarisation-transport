using LinearAlgebra, DifferentialEquations, Plots

# Constants
M = 1.0  # Mass of the black hole
a = 0.9

observer = (1000.0, deg2rad(90), 0.0, 0.0)

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

pl_Econserved = plot(
    xlabel="Affine Parameter λ",
    ylabel="|ΔE| (log scale)",
    title="Energy Difference |E₀ - E(λ)|",
    legend=:outerright,
    yscale=:log10,
)

for i in 1:length(b_values)
    λ_vals = λ_vals_all[i]
    abs_E_diff = abs.(E_vals_all[i] .- E_vals_all[i][1])
    abs_E_diff[abs_E_diff .== 0] .= 1e-10  # Replace zeros with small values

    plot!(
        pl_Econserved,
        λ_vals, abs_E_diff,
        label="|ΔE| for b=$(b_values[i])",
        lw=1.5
    )
end

display(pl_Econserved)


pl_Lconserved = plot(
    xlabel="Affine Parameter λ",
    ylabel="|ΔL| (log scale)",
    title="Momentum Difference |L₀ - L(λ)|",
    legend=:outerright,
    yscale=:log10,
)

for i in 1:length(b_values)
    λ_vals = λ_vals_all[i]
    abs_L_diff = abs.(L_vals_all[i] .- L_vals_all[i][1])
    abs_L_diff[abs_L_diff .== 0] .= 1e-10  # Replace zeros with small values

    plot!(
        pl_Lconserved,
        λ_vals, abs_L_diff,
        label="|ΔL| for b=$(b_values[i])",
        lw=1.5
    )
end

display(pl_Lconserved)

pl_Qconserved = plot(
    xlabel="Affine Parameter λ",
    ylabel="|ΔQ| (log scale)",
    title="Carter Const. Difference |Q₀ - Q(λ)|",
    legend=:outerright,
    yscale=:log10,
)

for i in 1:length(b_values)
    λ_vals = λ_vals_all[i]
    abs_Q_diff = abs.(Q_vals_all[i] .- Q_vals_all[i][1])
    abs_Q_diff[abs_Q_diff .== 0] .= 1e-10  # Replace zeros with small values

    plot!(
        pl_Qconserved,
        λ_vals, abs_Q_diff,
        label="|ΔQ| for b=$(b_values[i])",
        lw=1.5
    )
end

display(pl_Qconserved)