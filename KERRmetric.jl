using DifferentialEquations
using Plots

# Kerr metric parameters (constants)
M = 1.0    # Mass of the black hole
a = 0.5    # Spin parameter of the black hole
E = 1.0    # Energy of the photon
L_z = 2.0  # Angular momentum about z-axis
Q = 1.0    # Carter constant

# Define the Δ and Σ functions in the Kerr metric
function Δ(r)
    return r^2 - 2 * M * r + a^2
end

# Radial part of the geodesic equation, R(r)
function R(r)
    term1 = E * (r^2 + a^2) - a * L_z
    term2 = r^2 + (L_z - a * E)^2 + Q
    return term1^2 - Δ(r) * term2
end

# Angular part of the geodesic equation, Θ(θ)
function Θ(θ)
    term1 = a^2 * (1 - E^2)
    term2 = L_z^2 / sin(θ)^2
    return Q - cos(θ)^2 * (term1 + term2)
end

# Derivative of the radial part with respect to r, dR/dr
function dR_dr(r)
    term1 = 2 * (E * (r^2 + a^2) - a * L_z) * E * 2 * r
    term2 = - (2 * r - 2 * M) * (r^2 + (L_z - a * E)^2 + Q)
    term3 = Δ(r) * (2 * r)
    return term1 + term2 - term3
end

# Derivative of the angular part with respect to θ, dΘ/dθ
function dΘ_dθ(θ)
    term1 = sin(θ) * cos(θ) * (a^2 * (1 - E^2) + L_z^2 / sin(θ)^2)
    term2 = cos(θ)^2 * (L_z^2 * cos(θ) / sin(θ)^3)
    return -2 * term1 - 2 * term2
end

# Define the system of ODEs for geodesic motion
function geodesic_equations!(du, u, params, λ)
    r, θ, φ, t, pr, pθ = u  # unpack state variables

    # Compute Δ for the given r
    Δ_r = Δ(r)

    # Equations for dr/dλ and dθ/dλ
    du[1] = pr  # dr/dλ
    du[2] = pθ  # dθ/dλ
    du[3] = (L_z / sin(θ)^2) - a * E + a * (E * (r^2 + a^2) - a * L_z) / Δ_r  # dφ/dλ
    du[4] = E * (r^2 + a^2) / Δ_r - a * L_z / Δ_r  # dt/dλ

    # Radial and angular equations of motion
    du[5] = -0.5 * dR_dr(r)  # d(pr)/dλ
    du[6] = -0.5 * dΘ_dθ(θ)  # d(pθ)/dλ
end

# Initial conditions
r0 = 10.0
θ0 = π / 2
φ0 = 0.0
t0 = 0.0

# Initial momentum components (based on R and Θ)
pr0 = sqrt(R(r0))    # dr/dλ
pθ0 = sqrt(Θ(θ0))    # dθ/dλ

# Initial state vector
u0 = [r0, θ0, φ0, t0, pr0, pθ0]

# Time span (λ parameter)
λ_span = (0.0, 50.0)

# Solve the ODE
prob = ODEProblem(geodesic_equations!, u0, λ_span)
sol = solve(prob, Tsit5())

# Plot results
plot(sol, vars=(1, 2), xlabel="r", ylabel="θ", title="Photon Path in Kerr Metric")
