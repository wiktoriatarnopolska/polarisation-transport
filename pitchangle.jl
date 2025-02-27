using LinearAlgebra
using DifferentialEquations
using Plots

g = metric(r_hit, θ_hit)

g_tt = g[1,1]
g_tϕ = g[1,4]
g_rr = g[2,2]
g_θθ = g[3,3]
g_ϕϕ = g[4,4]

inv_g = inv(g)

gtt = inv_g[1,1]
gtϕ = inv_g[1,4]
grr = inv_g[2,2]
gθθ = inv_g[3,3]
gϕϕ = inv_g[4,4]

E, L, _ = calculate_conserved_quantities(g, kμ, a, θ_hit)
μ = 84.99398712522732

μ_rad = deg2rad(μ)

e = lnrf_tetrad(r_hit, θ_hit, a)

pϕ_LNRF = e[4,4] * (-gtϕ * E + gϕϕ * L) + e[4, 1] * (-gtϕ * E + gϕϕ * L)

Q = pϕ_LNRF ^2 * (cos(μ_rad))^2

W = 2 * gtϕ * E * L - gtt * E^2 - gϕϕ * L^2

p_θ_ = sqrt((Q + (cos(μ_rad))^2 * W * grr^2) / (gθθ * (gθθ^2 * (sin(μ_rad))^2 + (cos(μ_rad))^2 * grr^2)))
p_r_ = sqrt((W - gθθ * p_θ_^2)/(grr))

# Compute momenta magnitudes (absolute values)
p_θ_ = sqrt((Q + (cos(μ_rad))^2 * W * grr^2) / (gθθ * (gθθ^2 * (sin(μ_rad))^2 + (cos(μ_rad))^2 * grr^2)))
p_r_ = sqrt((W - gθθ * p_θ_^2) / grr)

# Determine the correct signs based on movement towards the observer
p_r_sign = sign(r0 - r_hit)  # Positive if moving outward, negative if moving inward
p_θ_sign = sign(θ0 - θ_hit)  # Positive if moving upwards, negative if moving downwards

# fixing the sign + raising the index
p_r_upp = grr * p_r_sign * p_r_
p_θ_upp = gθθ * p_θ_sign * p_θ_

p_t_upp = -gtt * E + gtϕ * L
p_ϕ_upp = -gtϕ * E + gϕϕ * L

p0_rev = [λ0, r_hit, θ_hit, ϕ_hit, p_t_upp, p_r_upp, p_θ_upp, p_ϕ_upp]

tspan_rev = (0.0, 5000.0)
# Create a ContinuousCallback that uses rootfinding (which by default uses a bisection‐like procedure)
obs_callback = ContinuousCallback(observer_condition, observer_affect!;
                                  rootfind = true,
                                  abstol = 1e-9, reltol = 1e-9)

prob_rev = ODEProblem(intprob!, p0_rev, tspan_rev)
sol_rev = solve(prob_rev, Tsit5(), callback = obs_callback, abstol=1e-9, reltol=1e-9)

λ_vals = []
E_vals = []
L_vals = []
Q_vals = []
H_vals = []

# Check conservation of quantities at each time step
for i in 1:length(sol_rev)
    # Extract position and momentum from the solution
    r = sol_rev[i][2]
    θ = sol_rev[i][3]
    ϕ = sol_rev[i][4]
    p_t = sol_rev[i][5]
    p_r = sol_rev[i][6]
    p_θ = sol_rev[i][7]
    p_ϕ = sol_rev[i][8]
    p = [p_t, p_r, p_θ, p_ϕ]

    # Metric at the current position
    g = metric(r, θ)

    
    # Calculate conserved quantities at current step
    E, L, Q = calculate_conserved_quantities(g, p, a, θ)
    H = g[1,1]*p_t^2 + g[2,2]*p_r^2 + g[3,3]*p_θ^2 + g[4,4]*p_ϕ^2 + 2*g[1,4]*p_t*p_ϕ


    # Store the values for plotting
    push!(E_vals, E)
    push!(L_vals, L)
    push!(Q_vals, Q)
    push!(H_vals, H)
end

# Extract positions and times
t_vals_rev = sol_rev.t
r_vals_rev = [sol_rev[i][2] for i in 1:length(sol_rev)]
θ_vals_rev = [sol_rev[i][3] for i in 1:length(sol_rev)]
ϕ_vals_rev = [sol_rev[i][4] for i in 1:length(sol_rev)]

x_vals_rev = [sqrt(r^2 + a^2) * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals_rev, θ_vals_rev, ϕ_vals_rev)]
y_vals_rev = [sqrt(r^2 + a^2) * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals_rev, θ_vals_rev, ϕ_vals_rev)]
z_vals_rev = [r * cos(θ) for (r, θ) in zip(r_vals_rev, θ_vals_rev)]


disc_hits_plot(r_in, r_out, r_horizon, x_hits, y_hits, x_vals, y_vals)
plot!(
    x_vals_rev, y_vals_rev,
    lw = 1.0,
    linecolor = :red,
    label = "returning ray",
    linestyle=:dash
)

conservations_plot(sol_rev, E_vals, L_vals, Q_vals, H_vals, false)
