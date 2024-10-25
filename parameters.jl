# Parameters for Kerr black hole

M = 10.0 * sol_M  # Mass of the black hole (kg)
f_col = 1.7  # Hardening factor

# redundant
#gravitational_radius(r_g::Float64) = r_g
#r_g = G * M / c^2  # Gravitational radius (meters)

# Observer inclination
θobs = deg2rad(30)

export M, am, f_col, r_g, rin, rout, θobs