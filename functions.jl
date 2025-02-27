# Function to calculate ISCO for a Kerr black hole (Bardeen et al. 1972)
function isco_radius(am::Float64)
    Z1 = 1 + (1 - am^2)^(1/3) * ((1 + am)^(1/3) + (1 - am)^(1/3))
    Z2 = sqrt(3 * am^2 + Z1^2)

    if am >= 0
        return (3 + Z2 - sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)))  # In units of r_g
    else
        return (3 + Z2 + sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)))  # In units of r_g
    end
end

# Function to calculate Novikov-Thorne radial profile (Page+Thorne 1974)
function novikov_thorne_profile(r::Float64, isco::Float64,M::Float64, am::Float64, M_dot::Float64, f_col::Float64)
    ξ_1 = 2 * cos((acos(am))/3 - π/3)
    ξ_2 = 2 * cos((acos(am))/3 + π/3)
    ξ_3 = -2 * cos((acos(am))/3)
    ξ_ms = sqrt(isco)  # ISCO is already in units of r_g
    ξ = sqrt(r)

    # Correct formula for temperature profile
    f_ξa = (ξ^4 * (ξ^3 - 3 * ξ + 2 * am))^(-1) * 
    (
        ξ - ξ_ms - (3/2) * am * log(ξ / ξ_ms)
        - (3 * (ξ_1 - am)^2 * (ξ_1 * (ξ_1 - ξ_2) * (ξ_1 - ξ_3))^(-1) * log((ξ - ξ_1) / (ξ_ms - ξ_1)))
        - (3 * (ξ_2 - am)^2 * (ξ_2 * (ξ_2 - ξ_1) * (ξ_2 - ξ_3))^(-1) * log((ξ - ξ_2) / (ξ_ms - ξ_2)))
        - (3 * (ξ_3 - am)^2 * (ξ_3 * (ξ_3 - ξ_1) * (ξ_3 - ξ_2))^(-1) * log((ξ - ξ_3) / (ξ_ms - ξ_3)))
    )

    T_amr = 741 * f_col * (M / sol_M)^(-1/2) * (M_dot)^(1/4) * (f_ξa)^(1/4) # keV
    #T_K = T_amr * 1e3 / k_B
    
    return T_amr
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

# Kerr Metric components as a function of r and θ
function metric(r, θ)
    Δ = r^2 - 2 * r + a^2
    Σ = r^2 + a^2 * cos(θ)^2

    g_tt = -(1 - (2 * r/ Σ))
    g_rr = Σ/Δ
    g_θθ = Σ
    g_ϕϕ = (sin(θ)^2 / Σ) * ((r^2 + a^2)^2 - Δ * a^2 * sin(θ)^2)
    g_tϕ = -(2 * r * a * sin(θ)^2) / Σ

    g = [
        g_tt     0        0        g_tϕ;
        0        g_rr     0        0;
        0        0        g_θθ     0;
        g_tϕ     0        0        g_ϕϕ
    ]
    return g
end

function horizon(a)
    r_horizon = 1 + sqrt(1 - a^2)
    return r_horizon
end

function disc_hit_condition(u, t, integrator)
    r = u[2]
    θ = u[3]
    return abs(θ - π/2) < 1e-6 ? 0.0 : θ - π/2
    #return abs(θ - π/2)
end

function disc_hit_affect!(integrator)
    r = integrator.u[2]
    θ = integrator.u[3]
    if r_in <= r <= r_out
        #println("Hit detected: r = $r, θ = $θ")
        terminate!(integrator)
    end
end

# Transform the photon's (x, y, 0) coordinates to BH Cartesian coordinates
function transform_to_bh_coords(x, y, observer, a)
    r_obs, θ_obs, ϕ_obs = observer[1:3]

    D = (sqrt(r_obs^2 + a^2)) * sin(θ_obs) - y * cos(θ_obs)
    
    # Cartesian coordinates relative to the BH
    x_bh = D * cos(ϕ_obs) - x * sin(ϕ_obs)
    y_bh = D * sin(ϕ_obs) + x * cos(ϕ_obs)
    z_bh = r_obs * cos(θ_obs) + y * sin(θ_obs)
    
    return x_bh, y_bh, z_bh
end

# Convert BH Cartesian coordinates to Boyer-Lindquist coordinates
function to_boyer_lindquist(x_bh, y_bh, z_bh, a)

    δ = x_bh^2 + y_bh^2 + z_bh^2 - a^2

    r = sqrt((δ + sqrt(δ^2 + 4 * a^2 * z_bh^2))/2)
    θ = acos(z_bh / r)
    ϕ = atan(y_bh, x_bh)  # Note: atan2(y_bh, x_bh) can be used for a more robust implementation
    return r, θ, ϕ
end

# Solve for p^t given the metric and momentum components
function solve_pt(g, p_r, p_θ, p_ϕ)
    g_tt = g[1, 1]
    g_tϕ = g[1, 4]
    g_rr = g[2, 2]
    g_θθ = g[3, 3]
    g_ϕϕ = g[4, 4]

    # Quadratic coefficients for p^t
    A = g_tt
    B = g_tϕ * p_ϕ
    C = g_rr * p_r^2 + g_θθ * p_θ^2 + g_ϕϕ * p_ϕ^2

    Δ = (2 * B)^2 - 4 * A * C
    if Δ < 0
        error("No real solution for p^t. Check initial conditions.")
    end

    # Use the negative root for future-directed motion
    p_t = (-2 * B - sqrt(Δ)) / (2 * A)
    
    return p_t
end

# Function to calculate E and Lz from p^t, p^phi and metric components
function calculate_energy_angular_momentum(g, p_t, p_ϕ)
    g_tt = g[1, 1]
    g_tϕ = g[1, 4]
    g_ϕϕ = g[4, 4]

    # Calculate Energy E
    E = - (g_tt * p_t + g_tϕ * p_ϕ)

    # Calculate Angular Momentum L_z
    L_z = g_tϕ * p_t + g_ϕϕ * p_ϕ

    return E, L_z
end


# Function to calculate energy (E), angular momentum (Lz), and Carter constant (Q)
function calculate_conserved_quantities(g, p, a, θ)
    g_tt, g_tϕ, g_rr, g_θθ, g_ϕϕ = g[1, 1], g[1, 4], g[2, 2], g[3, 3], g[4, 4]

    # Extract momentum components
    #p_t, p_r, p_θ, p_ϕ = p

    # Calculate Energy E
    # E = - (g_tt * p_t + g_tϕ * p_ϕ)

    # # Calculate Angular Momentum Lz
    # L_z = g_tϕ * p_t + g_ϕϕ * p_ϕ

    # # Calculate Carter Constant Q
    # Q = (g_θθ * p_θ)^2 + cos(θ)^2 * ( - a^2 * (g_tt * p_t + g_tϕ * p_ϕ)^2 + (g_tϕ * p_t + g_ϕϕ * p_ϕ)^2 / sin(θ)^2)
    
    E = -g_tt * p[1] - g_tϕ * p[4] # Energy-like quantity
    L_z = g_tϕ * p[1] + g_ϕϕ * p[4]  # Angular momentum-like quantity
    Q = (g_θθ * p[3])^2 + cos(θ)^2 * ( - a^2 * (g_tt * p[1]+ g_tϕ * p[4])^2 + ((g_ϕϕ * p[4] + g_tϕ * p[1])^2 / sin(θ)^2))


    return E, L_z, Q

end

# Define the ODE function using an in-place form
function intprob!(du, u, p, λ)
    # Current position and momentum
    x = u[1:4]
    p = u[5:8]

    r_val, θ_val = x[2], x[3]

    # Initialize du
    du .= 0.0

    # Terminate the integration if r crosses the event horizon
    if r_val <= r_horizon
        du[1:4] .= p
        return
    end

    # Compute analytical Christoffel symbols
    Γ = compute_christoffel_analytical(r_val, θ_val)

    # Compute the derivatives of velocity
    dp = zeros(4)
    for μ in 1:4
        sum_ = 0.0
        for ν in 1:4
            for λ in 1:4
                sum_ += Γ[μ, ν, λ] * p[ν] * p[λ]
            end
        end
        dp[μ] = -sum_
    end

    # Assign derivatives to du
    du[1:4] .= p
    du[5:8] .= dp
end

# Condition function: it returns the difference between the current radial coordinate and the observer's radius.
function observer_condition(u, t, integrator)
    return u[2] - observer[1]
end

# Affect function: simply terminate the integration when the event is found.
function observer_affect!(integrator)
    terminate!(integrator)
end

function lnrf_tetrad(r, θ, a)
    # Calculate metric-related quantities
    Σ = r^2 + (a * cos(θ))^2
    Δ = r^2 - 2*r + a^2
    A = (r^2 + a^2)^2 - a^2 * Δ * sin(θ)^2
    ω = 2*a*r/A  # Frame dragging angular velocity

    # Initialize tetrad matrix (each row is a tetrad vector)
    e = zeros(4, 4)
    
    # Time-like tetrad vector (e^μ_0)
    e[1,1] = sqrt(A / (Δ * Σ))                                   # t component
    e[1,2] = 0                                                 # r component
    e[1,3] = 0                                                 # θ component
    e[1,4] = 2 * a * r/(sqrt(Δ*Σ*A))                             # φ component
    
    # Radial tetrad vector (e^μ_1)
    e[2,1] = 0
    e[2,2] = sqrt(Δ/Σ)
    e[2,3] = 0
    e[2,4] = 0
    
    # Theta tetrad vector (e^μ_2)
    e[3,1] = 0
    e[3,2] = 0
    e[3,3] = 1 / sqrt(Σ)
    e[3,4] = 0
    
    # Phi tetrad vector (e^μ_3)
    e[4,1] = 0
    e[4,2] = 0
    e[4,3] = 0
    e[4,4] = sqrt(Σ / A) / sin(θ)
    
    return e
end

function incidence_angle(r_hit, θ_hit, ϕ_hit, p)
    
    # Get LNRF tetrad at hit point
    e = lnrf_tetrad(r_hit, θ_hit, a)
    
    # Transform momentum to LNRF frame
    p_lnrf = zeros(4)
    for μ in 1:4
        for ν in 1:4
            p_lnrf[μ] += e[μ,ν] * p[ν]
        end
    end
    
    # In the LNRF frame, the disc's normal vector at θ = π/2 is [0, 0, 1, 0]
    # Only need the spatial components for the angle calculation
    p_spatial = p_lnrf[2:4]
    n_disc = [0.0, 1.0, 0.0]  # Normal vector in LNRF frame
    
    norm_p_spatial = sqrt(sum(p_spatial .^ 2))
    norm_n_disc = sqrt(sum(n_disc .^ 2))

    # Calculate angle using dot product
    cos_angle = abs(dot(p_spatial, n_disc)) / (norm_p_spatial * norm_n_disc)
    angle = acos(cos_angle)
    
    return rad2deg(angle)  # Convert to degrees
end

export isco_radius, 
novikov_thorne_profile, 
compute_christoffel_analytical, 
metric, 
horizon, 
disc_hit_condition, 
disc_hit_affect!,
transform_to_bh_coords,
to_boyer_lindquist,
solve_pt,
calculate_energy_angular_momentum,
calculate_conserved_quantities,
intprob!,
observer_condition,
observer_affect!,
lnrf_tetrad,
incidence_angle