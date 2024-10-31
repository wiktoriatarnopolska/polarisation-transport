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
        - (3 * (ξ_1 - am)^2 * (ξ_1 * (ξ_1 - ξ_2) * (ξ_1 - ξ_3))^(-1) * log((ξ - ξ_1) / (ξ_ms - ξ_2)))
        - (3 * (ξ_3 - am)^2 * (ξ_3 * (ξ_3 - ξ_1) * (ξ_3 - ξ_2))^(-1) * log((ξ - ξ_3) / (ξ_ms - ξ_3)))
    )

    T_amr = 741 * f_col * (M / sol_M)^(-1/2) * (M_dot / sol_M)^(1/4) * (f_ξa)^(1/4) # keV
    T_K = T_amr * 1e3 / k_B
    
    return T_K
end

export isco_radius, novikov_thorne_profile 