using Roots
using Plots

# Constants
gravitational_radius(r_g::Float64) = r_g
solar_mass(sol_M::Float64) = sol_M
accretion_rate(M_dot::Float64) = M_dot
gravitational_constant(G::Float64) = G
speed_of_light(c::Float64) = c

G = 6.6743e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 2.99792458e8  # Speed of light (m/s)
sol_M = 1.989e+30  # Solar mass (kg)
M_dot = 1.4e18  # Accretion rate (kg/s)

# Parameters for Kerr black hole
M = 3.0 * sol_M  # Mass of the black hole (kg)
r_g = G * M / c^2  # Gravitational radius (meters)
am = 0.998  # Angular momentum a/M
f_col = 1.7  # Hardening factor

#######################################################################################################

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

# Calculate ISCO (in units of r_g)
isco = isco_radius(am)
println("ISCO radius: $isco r_g")

# Function to calculate Novikov-Thorne radial profile (Page+Thorne 1974)
function novikov_thorne_profile(ξ::Float64, M::Float64, am::Float64, M_dot::Float64, f_col::Float64)
    ξ_1 = 2 * cos(acos(am) - π/3)
    ξ_2 = 2 * cos(acos(am) + π/3)
    ξ_3 = -2 * cos(acos(am))
    ξ_ms = sqrt(isco)  # ISCO is already in units of r_g

    # Correct formula for temperature profile
    f_ξa = (ξ^4 * (ξ^3 - 2 * ξ + 2 * am))^(-1) * (
        ξ - ξ_ms - (3/2) * am * log(ξ / ξ_ms)
        - (3 * (ξ_1 - am)^2 * (ξ_1 * (ξ_1 - ξ_2) * (ξ_1 - ξ_3))^(-1) * log((ξ - ξ_1) / (ξ_ms - ξ_2)))
        - (3 * (ξ_3 - am)^2 * (ξ_3 * (ξ_3 - ξ_1) * (ξ_3 - ξ_2))^(-1) * log((ξ - ξ_3) / (ξ_ms - ξ_3)))
    )

    T_amr = 741 * f_col * (M / sol_M)^(-1/2) * (M_dot / sol_M)^(1/4) * (f_ξa)^(1/4)
    
    return T_amr
end

# Generate radii from ISCO to 1000 r_g (radii are now dimensionless, in units of r_g)
radii = range(isco, stop=1000, length=100)  # Radii in units of r_g

# Calculate Novikov-Thorne radial profile for each radius
profiles = [novikov_thorne_profile(r, M, am, M_dot, f_col) for r in radii]

# Plotting the radial profile
plot(radii, profiles, xlabel="Radius (r_g)", ylabel="Temperature (K)", title="Novikov-Thorne Radial Profile", legend=false)

# Marking the ISCO on the plot
vline!([isco], label="ISCO", linestyle=:dash, color=:red)

