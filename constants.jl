# Constants

solar_mass(sol_M::Float64) = sol_M
accretion_rate(M_dot::Float64) = M_dot
gravitational_constant(G::Float64) = G
speed_of_light(c::Float64) = c
Boltzmann_constant(k_B::Float64) = k_B
Planck_constant(h::Float64) = h

# GR regime
G = 1
c = 1
sol_M = 1.989e+30  # Solar mass (kg)
M_dot = 10e19  # Accretion rate (kg/s)
h = 4.135667696e-15  # Planck constant (eV·s)
k_B = 8.617333262145e-5  # Boltzmann constant (eV/K)

export G, c, sol_M, M_dot, h, k_B 