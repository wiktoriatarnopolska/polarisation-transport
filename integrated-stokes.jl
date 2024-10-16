# Constants and Parameters
#const ISCO = 6.0  # ISCO in terms of r_g for non-spinning Kerr
const ISCO = 1.0   # ISCO in terms of r_g for spinning Kerr BH
const r_max = 1000.0  # Outer radius of the disk in terms of r_g
const f_col = 1.8  # Hardening factor
const M = 14.0  # Mass of the black hole in solar masses
const G = 6.67430e-11  # Gravitational constant in SI units
const c = 2.998e8  # Speed of light in m/s
const sigma = 5.670374419e-8  # Stefan-Boltzmann constant in SI

# Functions to calculate temperature profile, flux, Chandrasekhar profile, Stokes parameters

# Gravitational radius (r_g) in meters
function r_g(M_solar)
    return G * M_solar/ c^2
end

# Novikov-Thorne temperature profile (assuming steady state)
function temperature_profile(r, M_dot)
    r_g_m = r_g(M)  # Gravitational radius in meters
    r_phys = r * r_g_m  # Physical radius
    T = f_col * (M_dot * G * M * 1.989e30 / (8 * π * sigma * r_phys^3))^0.25  # Temperature in Kelvin
    return T
end

# Calculate the flux at a given radius based on temperature
function flux_at_radius(T)
    return sigma * T^4  # Flux in W/m^2
end

# Chandrasekhar limb-darkening profile at optical depth = infinity
function chandrasekhar_profile(mu)
    return 1.5 * mu  # Simple limb-darkening formula
end

# Define local Stokes parameters (I, Q, U, V) based on flux
function stokes_parameters(F, angle)
    I = F  # Intensity proportional to flux
    Q = F * cos(2 * angle)
    U = F * sin(2 * angle)
    V = 0  # Assuming no circular polarization in this simple model
    return I, Q, U, V
end

# Polarization angle based on local disk properties (simplified here)
function polarization_angle(r)
    return atan(r / r_max)  # Simplified angle formula
end

# Integrate Stokes parameters over the disk
function integrate_stokes(r_range, M_dot)
    I_total, Q_total, U_total, V_total = 0.0, 0.0, 0.0, 0.0
    
    for r in r_range
        T = temperature_profile(r, M_dot)
        F = flux_at_radius(T)
        angle = polarization_angle(r)
        I, Q, U, V = stokes_parameters(F, angle)
        
        # Summing over the disk
        I_total += I
        Q_total += Q
        U_total += U
        V_total += V
    end
    
    return I_total, Q_total, U_total, V_total
end

# Main program execution

# Range of radii from ISCO to r_max
r_range = range(ISCO, stop=r_max, length=1000)

# Mass accretion rate (in kg/s, approximate for a supermassive BH)
M_dot = 1.0e18  # Example accretion rate

# Integrate the Stokes parameters over the disk
I_total, Q_total, U_total, V_total = integrate_stokes(r_range, M_dot)

# Print final results
println("Integrated Stokes Parameters:")
println("I_total: $I_total, Q_total: $Q_total, U_total: $U_total, V_total: $V_total")
