# Parameters for Kerr black hole

M = 10.0 * sol_M  # Mass of the black hole (kg)
f_col = 1.8  # Hardening factor

# redundant
#gravitational_radius(r_g::Float64) = r_g
#r_g = G * M / c^2  # Gravitational radius (meters)

export M, am, f_col, r_g, rin, rout, θobs