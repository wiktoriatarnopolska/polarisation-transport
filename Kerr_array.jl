using LinearAlgebra, DifferentialEquations, Plots

# Constants
M = 1.0  # Mass of the black hole
a = 0.9

r_horizon = horizon(a)

observer = (1000.0, deg2rad(90), 0.0, 0.0)

tspan = (0.0, 5000.0)

# Define a range of impact parameters
# b_values = collect(-10.0:1:10.0)
x_values = collect(-10.0:1.0:10.0)  # x values from -10 to 10 with step size 1
x_values
deleteat!(x_values, 11)

# Arrays to store trajectories
x_vals_all = []
y_vals_all = []

λ_vals_all = []
E_vals_all = []
L_vals_all = []
Q_vals_all = []

# for b in b_values
for x in x_values
    # Compute initial position in Cartesian coordinates
    r0 = observer[1]
    θ0 = observer[2]
    ϕ0 = observer[3]
    λ0 = observer[4]

    # # Initial velocities in spherical coordinates
    # v_r = -1.0                 # Photon moving inward
    # v_θ = 0.0                  # Equatorial plane
    # v_ϕ = b / (r0^2)           # Angular velocity from impact parameter.
    #                            # this ensures that the photon starts with 
    #                            # an angular momentum corresponding to the 
    #                            # impact parameter b at the initial large radius.

    # # Compute the metric at the initial position
    # g0 = metric(r0, θ0)
    # g_tt = g0[1,1]
    # g_tϕ = g0[1,4]
    # g_rr = g0[2,2]
    # g_θθ = g0[3,3]
    # g_ϕϕ = g0[4,4]

    # # Solve quadratic equation for v_t
    # A = g_tt
    # B = g_tϕ * v_ϕ
    # C = g_rr * v_r^2 + g_θθ * v_θ^2 + g_ϕϕ * v_ϕ^2

    # Δ = (2 * B)^2 - 4 * A * C
    # if Δ < 0
    #     error("No real solution for v_t. Check initial conditions.")
    # end

    # v_t = (-2 * B - sqrt(Δ)) / (2 * A)  # Negative root for future-directed motion

    # # Complete the initial velocity vector
    # v = [v_t; v_r; v_θ; v_ϕ]

    y = 0

    # Transform to BH Cartesian coordinates
    x_bh, y_bh, z_bh = transform_to_bh_coords(x, y, observer, a)

    # Convert BH Cartesian coordinates to Boyer-Lindquist coordinates
    r, θ, ϕ = to_boyer_lindquist(x_bh, y_bh, z_bh, a)

    # Calculate initial radial, theta, and phi momentum components
    Σ = r^2 + a^2 * cos(θ)^2
    p_r = - (r * sqrt(r^2 + a^2) * sin(θ) * sin(observer[2]) * cos(ϕ - observer[3]) + (r^2 + a^2) * cos(θ) * cos(observer[2])) / Σ
    p_θ = (r * sin(θ) * cos(observer[2]) - sqrt(r^2 + a^2) * cos(θ) * sin(observer[2]) * cos(ϕ - observer[3])) / Σ
    p_ϕ = (sin(observer[2]) * sin(ϕ - observer[3])) / (sqrt(r^2 + a^2) * sin(θ))

    # Metric at the initial position
    g = metric(r, θ)

    # Solve for p^t using the null condition
    p_t = solve_pt(g, p_r, p_θ, p_ϕ)

    # Initial momentum vector
    p = [p_t, p_r, p_θ, p_ϕ]

    # Initial state vector
    u0 = [λ0, r0, θ0, ϕ0, p[1], p[2], p[3], p[4]]

    # Solve the ODE
    prob = ODEProblem(intprob!, u0, tspan)
    sol = solve(prob, Tsit5(), abstol=1e-9, reltol=1e-9)

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
for i in 1:length(x_values)
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
    yscale=:log10,
    legend=false
)

for i in 1:length(x_values)
    λ_vals = λ_vals_all[i]
    abs_E_diff = abs.(E_vals_all[i] .- E_vals_all[i][1])
    abs_E_diff[abs_E_diff .== 0] .= 1e-10  # Replace zeros with small values

    plot!(
        pl_Econserved,
        λ_vals, abs_E_diff,
        lw=1.5
    )
end

display(pl_Econserved)


pl_Lconserved = plot(
    xlabel="Affine Parameter λ",
    ylabel="|ΔL| (log scale)",
    title="Momentum Difference |L₀ - L(λ)|",
    yscale=:log10,
    legend=false
)

for i in 1:length(x_values)
    λ_vals = λ_vals_all[i]
    abs_L_diff = abs.(L_vals_all[i] .- L_vals_all[i][1])
    abs_L_diff[abs_L_diff .== 0] .= 1e-10  # Replace zeros with small values

    plot!(
        pl_Lconserved,
        λ_vals, abs_L_diff,
        lw=1.5
    )
end

display(pl_Lconserved)

pl_Qconserved = plot(
    xlabel="Affine Parameter λ",
    ylabel="|ΔQ| (log scale)",
    title="Carter Const. Difference |Q₀ - Q(λ)|",
    legend=:false,
    yscale=:log10,
)

for i in 1:length(x_values)
    λ_vals = λ_vals_all[i]
    abs_Q_diff = abs.(Q_vals_all[i] .- Q_vals_all[i][1])
    abs_Q_diff[abs_Q_diff .== 0] .= 1e-10  # Replace zeros with small values

    plot!(
        pl_Qconserved,
        λ_vals, abs_Q_diff,
        lw=1.5
    )
end

display(pl_Qconserved)

plt = plot(pl_Econserved, pl_Lconserved, pl_Qconserved, layout = (1, 3))
