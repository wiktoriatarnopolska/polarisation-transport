using LinearAlgebra
using DifferentialEquations
using Plots

# Constants
M = 1.0  # Mass of the black hole
a = 0.9

observer = (1000.0, deg2rad(60), 0.0, 0.0)

callback = ContinuousCallback(disc_hit_condition, disc_hit_affect!)

tspan = (0.0, 5000.0)

# Disc parameters
r_horizon = horizon(a)
r_in = isco_radius(a)       # Inner radius of the disc
r_out = 10.0                # Outer radius of the disc

# Arrays to store trajectories
x_vals_all = []
y_vals_all = []

# Initialize array to store intersection points
disc_hits = []

r0 = observer[1]
θ0 = observer[2]
ϕ0 = observer[3]
λ0 = observer[4]

# Transforming and solving for initial conditions
x, y = 6.0, 0.0 
x_bh, y_bh, z_bh = transform_to_bh_coords(x, y, observer, a)
r, θ, ϕ = to_boyer_lindquist(x_bh, y_bh, z_bh, a)

Σ = r^2 + a^2 * cos(θ)^2

# Initial radial, theta, phi momentum components
p_r = - (r * sqrt(r^2 + a^2) * sin(θ) * sin(θ0) * cos(ϕ - ϕ0) + (r^2 + a^2) * cos(θ) * cos(θ0))/ (Σ)
p_θ = (r * sin(θ) * cos(θ0) - sqrt(r^2 + a^2) * cos(θ) * sin(θ0) * cos(ϕ-ϕ0))/(Σ)
p_ϕ = (sin(θ0) * sin(ϕ-ϕ0))/(sqrt(r^2 + a^2) * sin(θ))

# Metric at the initial position
g = metric(r, θ)
g_tt = g[1,1]
g_tϕ = g[1,4]
g_rr = g[2,2]
g_θθ = g[3,3]
g_ϕϕ = g[4,4]

# Solve for p^t
p_t = solve_pt(g, p_r, p_θ, p_ϕ)

H = g[1,1]*p_t^2 + g[2,2]*p_r^2 + g[3,3]*p_θ^2 + g[4,4]*p_ϕ^2 + 2*g[1,4]*p_t*p_ϕ

# Initial momentum vector
p = [p_t, p_r, p_θ, p_ϕ]

# Update initial state vector for ODE solver
u0 = [λ0, r, θ, ϕ, p[1], p[2], p[3], p[4]]

# Calculate initial (?) Energy (E) and Angular Momentum (Lz)
E, L_z, Q = calculate_conserved_quantities(g, p, a, θ)

# Calculate the impact parameter b = Lz / E
b = L_z / E

# Continue with ODE solving using the updated u0
prob = ODEProblem(intprob!, u0, tspan)
sol = solve(prob, Tsit5(), callback = callback, abstol=1e-12, reltol=1e-12, dtmax=0.01)


# Extract positions and times
t_vals = sol.t
r_vals = [sol[i][2] for i in 1:length(sol)]
θ_vals = [sol[i][3] for i in 1:length(sol)]
ϕ_vals = [sol[i][4] for i in 1:length(sol)]

# Initialize flag and variables for disc intersection
hit_disc = false
r_hit = 0.0
ϕ_hit = 0.0
θ_hit = 0.0

p_t_hit = 0.0
p_r_hit = 0.0
p_θ_hit = 0.0
p_ϕ_hit = 0.0

# Check for intersection with the disc
for i in 1:length(θ_vals)-1
    θ1 = θ_vals[i]
    θ2 = θ_vals[i+1]

    # Check if θ crosses π/2 between θ1 and θ2
    if (θ1 - π/2) * (θ2 - π/2) <= 0
        # Interpolate to find the crossing point
        fraction = (π/2 - θ1) / (θ2 - θ1)
        r_cross = r_vals[i] + fraction * (r_vals[i+1] - r_vals[i])
        ϕ_cross = ϕ_vals[i] + fraction * (ϕ_vals[i+1] - ϕ_vals[i])
        θ_cross = θ1 + fraction * (θ2 - θ1)

        # Normalize ϕ_cross between 0 and 2π
        ϕ_cross = mod(ϕ_cross, 2π)

        #Interpolated momenta at the hit
        p_t_hit = sol[i][5] + fraction * (sol[i+1][5] - sol[i][5])
        p_r_hit = sol[i][6] + fraction * (sol[i+1][6] - sol[i][6])
        p_θ_hit = sol[i][7] + fraction * (sol[i+1][7] - sol[i][7])
        p_ϕ_hit = sol[i][8] + fraction * (sol[i+1][8] - sol[i][8])

        # λ_cross = sol.t[i] + fraction * (sol.t[i + 1] - sol.t[i])

        # p_t_hit = sol(λ_cross)[5]
        # p_r_hit = sol(λ_cross)[6]
        # p_θ_hit = sol(λ_cross)[7]
        # p_ϕ_hit = sol(λ_cross)[8]

        # Check if r_cross is within the disc's radial extent
        if r_in <= r_cross <= r_out
            hit_disc = true
            r_hit = r_cross
            ϕ_hit = ϕ_cross
            θ_hit = θ_cross

            # p_t_hit = p_t
            # p_r_hit = p_r
            # p_θ_hit = p_θ
            # p_ϕ_hit = p_ϕ
            break  # Exit loop after first intersection
        end
    end
end

# Record the intersection if it occurs
# Record the intersection with momenta
if hit_disc
    push!(disc_hits, (b, r_hit, ϕ_hit, p_t_hit, p_r_hit, p_θ_hit, p_ϕ_hit))
end

x_vals = [sqrt(r^2 + a^2) * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
y_vals = [sqrt(r^2 + a^2) * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
z_vals = [r * cos(θ) for (r, θ) in zip(r_vals, θ_vals)]

# Append trajectories
push!(x_vals_all, x_vals)
push!(y_vals_all, y_vals)

println("Number of disc hits: ", length(disc_hits))

# Use disc_hits to analyze or visualize the geodesics that hit the disc
for (b, r_hit, ϕ_hit) in disc_hits
    println("Geodesic with impact parameter b = $b hits the disc at (r, ϕ, θ) = ($r_hit, $ϕ_hit, $θ_hit)") 
end

#### PLOTTING ################################################################################

# Convert the disc hits to Cartesian coordinates
x_hits = [r_hit * sin(π/2) * cos(ϕ_hit) for (_, r_hit, ϕ_hit) in disc_hits]
y_hits = [r_hit * sin(π/2) * sin(ϕ_hit) for (_, r_hit, ϕ_hit) in disc_hits]
z_hits = [r_hit * cos(π/2)]
# Plot the disc boundaries
θ_values = range(0, 2π, length=500)
x_inner = [r_in * cos(θ) for θ in θ_values]
y_inner = [r_in * sin(θ) for θ in θ_values]
x_outer = [r_out * cos(θ) for θ in θ_values]
y_outer = [r_out * sin(θ) for θ in θ_values]

# Initialize the plot
pl_disc = plot(
    x_outer, y_outer,
    seriestype = :shape,
    fillcolor = :mediumorchid2,
    linecolor = :rebeccapurple,
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    title = "Disc Hits Visualization",
    label = "Disc",
    xlim = (-20, 20),
    ylim = (-20, 20),
)

# Fill the inner disc area
plot!(
    x_inner, y_inner,
    seriestype = :shape,
    fillcolor = :white,
    linecolor = :rebeccapurple,
    label = false
)

# Overlay the event horizon
circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]
plot!(
    circle_x, circle_y,
    lw = 2,
    color = :purple4,
    label = "Event horizon"
)

# Overlay the ISCO
isco_x = [r_in * cos(θ) for θ in 0:0.01:2π]
isco_y = [r_in * sin(θ) for θ in 0:0.01:2π]
plot!(
    isco_x, isco_y,
    lw = 2,
    color = :mediumslateblue,
    label = "ISCO"
)

# Overlay the hit points
scatter!(
    pl_disc,
    x_hits, y_hits,
    color = :tan1,
    markerstrokecolor = :black,
    markersize = 6,
    label = "Disc Hits",
)

# Plot the trajectory
plot!(
    pl_disc,
    x_vals, y_vals,
    lw = 1.0,
    linecolor = :lightsteelblue,
    label = false
)

# Display the plot
display(pl_disc)



###########################################################################
#   TO TRACE BACK

# state at hit point

g = metric(r_hit, θ_hit)
g_tt = g[1,1]
g_tϕ = g[1,4]
g_rr = g[2,2]
g_θθ = g[3,3]
g_ϕϕ = g[4,4]

# Extract momentum components at hit
p_hit = [p_t_hit, p_r_hit, p_θ_hit, p_ϕ_hit]

# Compute conserved quantities at hit point
E_hit, Lz_hit, Q_hit = calculate_conserved_quantities(g, p_hit, a, θ_hit)

H_hit = g[1,1]*p_t_hit^2 + g[2,2]*p_r_hit^2 + g[3,3]*p_θ_hit^2 + g[4,4]*p_ϕ_hit^2 + 2*g[1,4]*p_t_hit*p_ϕ_hit


println("Initial Energy E0: ", E)
println("Energy at hit point E_hit: ", E_hit)
println("ΔE = ", abs(E - E_hit))

println("Initial Momentum L0: ", L_z)
println("Momentum at hit point Lz_hitt: ", Lz_hit)
println("ΔLz_hit = ", abs(L_z - Lz_hit))

println("Initial Carter const. Q0: ", Q)
println("Carter const. at hit point Q_hit: ", Q_hit)
println("ΔQ = ", abs(Q - Q_hit))

println("Initial H0: ", H)
println("H at hit point H_hit: ", H_hit)
println("ΔH = ", abs(H - H_hit))

# Reversing
#p_t = -p_t_hit
p_r = -p_r_hit
p_θ = π/2 - p_θ_hit
p_ϕ = p_ϕ_hit + π

# Δ = r_hit^2 - 2 * r_hit + a^2
# Σ = r_hit^2 + a^2 * (cos(θ_hit))^2
# A = (r_hit^2 + a^2)^2 - a^2 * Δ * sin(θ_hit)^2
# ω = 2 * a * r_hit / A

# et_t = sqrt(((r_hit^2 + a^2)^2 * Δ / A)/Σ)
# et_ϕ = - ω * et_t
# er_r = sqrt(Δ/Σ)
# eθ_θ = 1 / sqrt(Σ)
# eϕ_t = - ω * sqrt(Σ/A) * sin(θ_hit)
# eϕ_ϕ = sqrt(Σ/A) * sin(θ_hit)

# # # B-L -> LNRF
# # p_θ = eθ_θ * p_θ_hit
# # p_ϕ = et_ϕ * p_ϕ_hit + eϕ_ϕ * p_ϕ_hit

# # # Inverse

# # p_θ_rev = - p_θ
# # p_ϕ_rev = - p_ϕ

# # # Transform back 
# # p_θ = eθ_θ * p_θ_rev
# # p_ϕ = eϕ_ϕ * p_ϕ_rev - et_ϕ * p_ϕ_rev

# # correction_θ = p_θ_hit / p_θ
# # correction_ϕ = p_ϕ_hit / p_ϕ

# # To LNRF
# p_r_LNRF = er_r * p_r_hit
# p_θ_LNRF = eθ_θ * p_θ_hit
# p_ϕ_LNRF = eϕ_t * p_t_hit + eϕ_ϕ * p_ϕ_hit

# # Reverse in LNRF
# p_r = p_r_hit
# p_θ = ((p_θ_hit + π) % π)
# p_ϕ = (p_ϕ_hit + π) % (2*π)

# # Back to Boyer-Lindquist
# p_r = p_r_rev_LNRF / er_r
# p_θ = p_θ_rev_LNRF / eθ_θ
# p_ϕ = (p_ϕ_rev_LNRF - eϕ_t * p_t) / eϕ_ϕ

# p_r = p_r_hit
p_t = solve_pt(g, p_r, p_θ ,p_ϕ)

norm = g[1,1]*p_t^2 + g[2,2]*p_r^2 + g[3,3]*p_θ^2 + g[4,4]*p_ϕ^2 + 2*g[1,4]*p_t*p_ϕ

p = [p_t, p_r, p_θ, p_ϕ]

E, L_z, Q = calculate_conserved_quantities(g, p, a, θ_hit)

u0_rev = [λ0, r_hit, θ_hit, ϕ_hit, p_t, p_r, p_θ, p_ϕ]

tspan_rev = (0.0, -5000.0)

prob_rev = ODEProblem(intprob!, u0_rev, tspan_rev)
sol_rev = solve(prob_rev, Tsit5(), abstol=1e-14, reltol=1e-14, dtmax=1.0)

t_vals = sol_rev.t
r_vals = [sol_rev[i][2] for i in 1:length(sol_rev)]
θ_vals = [sol_rev[i][3] for i in 1:length(sol_rev)]
ϕ_vals = [sol_rev[i][4] for i in 1:length(sol_rev)]

x_vals = [sqrt(r^2 + a^2) * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
y_vals = [sqrt(r^2+ a ^2) * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]

push!(x_vals_all, x_vals)
push!(y_vals_all, y_vals)

θ_values = range(0, 2π, length=500)
x_inner = [r_in * cos(θ) for θ in θ_values]
y_inner = [r_in * sin(θ) for θ in θ_values]
x_outer = [r_out * cos(θ) for θ in θ_values]
y_outer = [r_out * sin(θ) for θ in θ_values]

pl_disc = plot(
    x_outer, y_outer,
    seriestype = :shape,
    fillcolor = :mediumorchid2,
    linecolor = :rebeccapurple,
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    title = "Disc Hits Visualization",
    label = "Disc",
    xlim = (-20, 20),
    ylim = (-20, 20),
)
plot!(
    x_inner, y_inner,
    seriestype = :shape,
    fillcolor = :white,
    linecolor = :rebeccapurple,
    label = false
)

circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]
plot!(
    circle_x, circle_y,
    lw = 2,
    color = :purple4,
    label = "Event horizon"
)

isco_x = [r_in * cos(θ) for θ in 0:0.01:2π]
isco_y = [r_in * sin(θ) for θ in 0:0.01:2π]
plot!(
    isco_x, isco_y,
    lw = 2,
    color = :mediumslateblue,
    label = "ISCO"
)

scatter!(
    pl_disc,
    x_hits, y_hits,
    color = :tan1,
    markerstrokecolor = :black,
    markersize = 6,
    label = "Disc Hits",
)

plot!(
    pl_disc,
    x_vals, y_vals,
    lw = 1.0,
    linecolor = :lightsteelblue,
    label = false
)

display(pl_disc)
