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

# Arrays to store λ (affine parameter) and conserved quantities for plotting for hitting geodesics
λ_vals = []
E_vals = []
L_vals = []
Q_vals = []
H_vals = []

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
sol = solve(prob, callback = callback, Tsit5(), abstol=1e-9, reltol=1e-9)

# Extract λ (affine parameter) and quantities for plotting
λ_vals = sol.t
E_vals = []
L_vals = []
Q_vals = []
H_vals = []


# Check conservation of quantities at each time step
for i in 1:length(sol)
    # Extract position and momentum from the solution
    r = sol[i][2]
    θ = sol[i][3]
    ϕ = sol[i][4]
    p_t = sol[i][5]
    p_r = sol[i][6]
    p_θ = sol[i][7]
    p_ϕ = sol[i][8]
    p = [p_t, p_r, p_θ, p_ϕ]

    # Metric at the current position
    g = metric(r, θ)

    
    # Calculate conserved quantities at current step
    E, Lz, Q = calculate_conserved_quantities(g, p, a, θ)
    H = g[1,1]*p_t^2 + g[2,2]*p_r^2 + g[3,3]*p_θ^2 + g[4,4]*p_ϕ^2 + 2*g[1,4]*p_t*p_ϕ


    # Store the values for plotting
    push!(E_vals, E)
    push!(L_vals, Lz)
    push!(Q_vals, Q)
    push!(H_vals, H)
end

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

        # #Interpolated momenta at the hit
        # p_t_hit = sol[i][5] + fraction * (sol[i+1][5] - sol[i][5])
        # p_r_hit = sol[i][6] + fraction * (sol[i+1][6] - sol[i][6])
        # p_θ_hit = sol[i][7] + fraction * (sol[i+1][7] - sol[i][7])
        # p_ϕ_hit = sol[i][8] + fraction * (sol[i+1][8] - sol[i][8])

        λ_cross = sol.t[i] + fraction * (sol.t[i + 1] - sol.t[i])

        p_t_hit = sol(λ_cross)[5]
        p_r_hit = sol(λ_cross)[6]
        p_θ_hit = sol(λ_cross)[7]
        p_ϕ_hit = sol(λ_cross)[8]

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

disc_hits_plot(r_in, r_out, r_horizon, x_hits, y_hits, x_vals, y_vals)

conservations_plot(sol, E_vals, L_vals, Q_vals, H_vals)



###########################################################################
#   TO TRACE BACK

# Arrays to store λ (affine parameter) and conserved quantities for plotting for hitting geodesics
#λ_vals = []
E_vals = []
L_vals = []
Q_vals = []
H_vals = []

# state at hit point

g = metric(r_hit, θ_hit)
g_tt = g[1,1]
g_tϕ = g[1,4]
g_rr = g[2,2]
g_θθ = g[3,3]
g_ϕϕ = g[4,4]

# Extract momentum components at hit
p_hit = [p_t_hit, p_r_hit, p_θ_hit, p_ϕ_hit]

# Reversing
p_t = p_t_hit
p_r = p_r_hit
p_θ = p_θ_hit
p_ϕ = p_ϕ_hit

H = g[1,1]*p_t^2 + g[2,2]*p_r^2 + g[3,3]*p_θ^2 + g[4,4]*p_ϕ^2 + 2*g[1,4]*p_t*p_ϕ

p = [p_t, p_r, p_θ, p_ϕ]

u0_rev = [λ0, r_hit, θ_hit, ϕ_hit, p_t, p_r, p_θ, p_ϕ]

# tspan_rev = (0.0, -λ_vals[end])

# Condition function: it returns the difference between the current radial coordinate and the observer's radius.
function observer_condition(u, t, integrator)
    return u[2] - observer[1]
end

# Affect function: simply terminate the integration when the event is found.
function observer_affect!(integrator)
    terminate!(integrator)
end

# Create a ContinuousCallback that uses rootfinding (which by default uses a bisection‐like procedure)
obs_callback = ContinuousCallback(observer_condition, observer_affect!;
                                  rootfind = true,
                                  abstol = 1e-9, reltol = 1e-9)

prob_rev = ODEProblem(intprob!, u0_rev, tspan_rev)
sol_rev = solve(prob_rev, Tsit5(), callback = obs_callback, abstol=1e-9, reltol=1e-9)

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

disc_hits_plot(r_in, r_out, r_horizon, x_hits, y_hits, x_vals, y_vals)
plot!(
    x_vals, y_vals,
    lw = 1.0,
    linecolor = :red,
    label = "returning ray",
    linestyle=:dash
)

conservations_plot(sol_rev, E_vals, L_vals, Q_vals, H_vals)
