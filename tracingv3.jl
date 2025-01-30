using LinearAlgebra
using DifferentialEquations
using Plots

# Constants
M = 1.0  # Mass of the black hole
a = 0.9

observer = (1000.0, deg2rad(60), 0.0, 0.0)

callback = ContinuousCallback(disc_hit_condition, disc_hit_affect!)

tspan = (0.0, 5000.0)

# Define a range of x values in the observer's image plane
#x_values = collect(-10.0:1.0:10.0)  # x values from -10 to 10 with step size 1
x_values = 6.0

# Disc parameters
r_horizon = horizon(a)
r_in = isco_radius(a)       # Inner radius of the disc
r_out = 10.0                  # Outer radius of the disc

# Arrays to store trajectories for photons that hit the disc
x_vals_hits = []
y_vals_hits = []

# Initialize array to store intersection points for photons that hit the disc
disc_hits = []

# Arrays to store λ (affine parameter) and conserved quantities for plotting for hitting geodesics
λ_vals_all = []
E_vals_all = []
L_vals_all = []
Q_vals_all = []

# Loop over all x values to compute photon trajectories
for x in x_values
    y = 0.0  # Could also vary y for more coverage

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

    # Calculate initial conserved quantities
    E_initial, Lz_initial, Q_initial = calculate_conserved_quantities(g, p, a)

    # Update initial state vector for ODE solver
    u0 = [observer[4], r, θ, ϕ, p[1], p[2], p[3], p[4]]

    # Define the ODE problem for the photon
    prob = ODEProblem(intprob!, u0, tspan)

    # Solve the ODE
    sol = solve(prob, Tsit5(), callback=callback, abstol=1e-12, reltol=1e-12, dtmax=1.0)

    # Extract λ (affine parameter) and quantities for plotting
    λ_vals = sol.t
    E_vals = []
    L_vals = []
    Q_vals = []

    # Check conservation of quantities at each time step
    for i in 1:length(sol)
        # Extract position and momentum from the solution
        r = sol[i][2]
        θ = sol[i][3]
        ϕ = sol[i][4]
        p = sol[i][5:8]

        # Metric at the current position
        g = metric(r, θ)

        # Calculate conserved quantities at current step
        E, Lz, Q = calculate_conserved_quantities(g, p, a)

        # Store the values for plotting
        push!(E_vals, E)
        push!(L_vals, Lz)
        push!(Q_vals, Q)
    end

    # Extract radial and azimuthal values for plotting
    r_vals = [sol[i][2] for i in 1:length(sol)]
    θ_vals = [sol[i][3] for i in 1:length(sol)]
    ϕ_vals = [sol[i][4] for i in 1:length(sol)]

    # Initialize flag and variables for disc intersection
    hit_disc = false
    r_hit = 0.0
    ϕ_hit = 0.0

    # Check for intersection with the disc
    for i in 1:length(θ_vals) - 1
        θ1 = θ_vals[i]
        θ2 = θ_vals[i + 1]

        # Check if θ crosses π/2 between θ1 and θ2
        if (θ1 - π/2) * (θ2 - π/2) <= 0
            # Interpolate to find the crossing point
            fraction = (π/2 - θ1) / (θ2 - θ1)
            r_cross = r_vals[i] + fraction * (r_vals[i + 1] - r_vals[i])
            ϕ_cross = ϕ_vals[i] + fraction * (ϕ_vals[i + 1] - ϕ_vals[i])

            # Normalize ϕ_cross between 0 and 2π
            ϕ_cross = mod(ϕ_cross, 2π)

            # Check if r_cross is within the disc's radial extent
            if r_in <= r_cross <= r_out
                hit_disc = true
                r_hit = r_cross
                ϕ_hit = ϕ_cross
                break  # Exit loop after first intersection
            end
        end
    end

    # Record the intersection and store the trajectory if it occurs
    if hit_disc
        # Store the (x, r_hit, ϕ_hit) for later use
        push!(disc_hits, (x, r_hit, ϕ_hit))

        # Append λ, E, L, Q for hitting geodesics only
        push!(λ_vals_all, λ_vals)
        push!(E_vals_all, E_vals)
        push!(L_vals_all, L_vals)
        push!(Q_vals_all, Q_vals)

        # Extract positions for plotting
        x_vals = [sqrt(r^2 + a^2) * sin(θ) * cos(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]
        y_vals = [sqrt(r^2 + a^2) * sin(θ) * sin(ϕ) for (r, θ, ϕ) in zip(r_vals, θ_vals, ϕ_vals)]


        # Append trajectories for photons that hit the disc
        push!(x_vals_hits, x_vals)
        push!(y_vals_hits, y_vals)
    end
end

# Print the number of disc hits
println("Number of disc hits: ", length(disc_hits))

# Use disc_hits to analyze or visualize the geodesics that hit the disc
for (x, r_hit, ϕ_hit) in disc_hits
    println("Geodesic with initial x = $x hits the disc at (r, ϕ) = ($r_hit, $ϕ_hit)")
end

# Plotting Conserved Quantities for Geodesics that Hit the Disc

# Plot Energy Difference |ΔE|
pl_Econserved = plot(
    xlabel = "Affine Parameter λ",
    ylabel = "|ΔE| (log scale)",
    title = "Energy Difference |E₀ - E(λ)|",
    legend = :outerright,
    yscale = :log10,
)

for i in 1:length(λ_vals_all)
    λ_vals = λ_vals_all[i]
    abs_E_diff = abs.(E_vals_all[i] .- E_vals_all[i][1])
    abs_E_diff[abs_E_diff .== 0] .= 1e-10  # Replace zeros with small values

    plot!(
        pl_Econserved,
        λ_vals, abs_E_diff,
        label = "|ΔE| for x=$(disc_hits[i][1])",
        lw = 1.5
    )
end

display(pl_Econserved)

# Plot Angular Momentum Difference |ΔL|
pl_Lconserved = plot(
    xlabel = "Affine Parameter λ",
    ylabel = "|ΔL| (log scale)",
    title = "Momentum Difference |L₀ - L(λ)|",
    legend = :outerright,
    yscale = :log10,
)

for i in 1:length(λ_vals_all)
    λ_vals = λ_vals_all[i]
    abs_L_diff = abs.(L_vals_all[i] .- L_vals_all[i][1])
    abs_L_diff[abs_L_diff .== 0] .= 1e-10  # Replace zeros with small values

    plot!(
        pl_Lconserved,
        λ_vals, abs_L_diff,
        label = "|ΔL| for x=$(disc_hits[i][1])",
        lw = 1.5
    )
end

display(pl_Lconserved)

# Plot Carter Constant Difference |ΔQ|
pl_Qconserved = plot(
    xlabel = "Affine Parameter λ",
    ylabel = "|ΔQ| (log scale)",
    title = "Carter Const. Difference |Q₀ - Q(λ)|",
    legend = :outerright,
    yscale = :log10,
)

for i in 1:length(λ_vals_all)
    λ_vals = λ_vals_all[i]
    abs_Q_diff = abs.(Q_vals_all[i] .- Q_vals_all[i][1])
    abs_Q_diff[abs_Q_diff .== 0] .= 1e-10  # Replace zeros with small values

    plot!(
        pl_Qconserved,
        λ_vals, abs_Q_diff,
        label = "|ΔQ| for x=$(disc_hits[i][1])",
        lw = 1.5
    )
end

display(pl_Qconserved)
