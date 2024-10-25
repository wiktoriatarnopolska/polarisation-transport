# Parameters for Kerr black hole

M = 3.0 * sol_M  # Mass of the black hole (kg)
am = 0.998 # Angular momentum a/M
f_col = 1.7  # Hardening factor

gravitational_radius(r_g::Float64) = r_g
r_g = G * M / c^2  # Gravitational radius (meters)

rin = isco_radius(am)
rout = 10

# Observer inclination
θobs = deg2rad(30)

export M, am, f_col, r_g, rin, rout, θobs