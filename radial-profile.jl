using Roots

# Constants
gravitational_radius(M::Float64) = M

# Function to calculate ISCO for a Kerr black hole
function isco_radius(a::Float64, M::Float64)
    if abs(a) > M
        error("Spin parameter |a| cannot be greater than the mass M.")
    end

    Z1 = 1 + (1 - (a/M)^2)^(1/3) * ((1 + a/M)^(1/3) + (1 - a/M)^(1/3))
    Z2 = sqrt(3 * (a/M)^2 + Z1^2)

    if a >= 0
        return M * (3 + Z2 - sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)))
    else
        return M * (3 + Z2 + sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)))
    end
end

# Function to calculate Novikov-Thorne radial profile
function novikov_thorne_profile(r::Float64, M::Float64, a::Float64, Mdot::Float64, f_col::Float64)
    # Constants
    r_g = gravitational_radius(M)
    
    # Relativistic functions for Novikov-Thorne model
    L_r = r^(-3/4)  # Placeholder for the actual function L(r, a)
    C_r = r^(-3/4)  # Placeholder for the actual function C(r, a)
    R_r = r^(-3/4) * (L_r / C_r)^(1/4)
    
    # Temperature profile based on Novikov-Thorne model
    T_r = 741 * f_col * (Mdot / (M * 1e18))^(1/4) * (M)^(-1/2) * R_r
    
    return T_r
end

# Parameters for Kerr black hole
M = 1.0            # Mass of the black hole (in units where G = c = 1)
a = 0.0 * M        # Spin parameter of the black hole (must be -M <= a <= M)
Mdot = 1.4e18      # Accretion rate in g/s
f_col = 1.7        # Hardening factor

# Calculate ISCO
isco = isco_radius(a, M)
println("ISCO radius: $isco r_g")

# Generate radii from ISCO to 1000 gravitational radii
radii = range(isco, stop=1000 * gravitational_radius(M), length=100)
println("Radii from ISCO to 1000 r_g:")
println(radii)

# Calculate Novikov-Thorne radial profile for each radius
profiles = [novikov_thorne_profile(r, M, a, Mdot, f_col) for r in radii]
println("Novikov-Thorne radial profile:")
println(profiles)

# Plot the radial profile
plot(radii, profiles, xlabel="Radius (r_g)", 
ylabel="Temperature (K)", title="Novikov-Thorne Radial Profile", 
legend=false)

vline!([isco], label="ISCO", linestyle=:dash, color=:red)