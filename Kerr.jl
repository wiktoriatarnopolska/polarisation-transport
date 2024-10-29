using LinearAlgebra, DifferentialEquations, Plots, ForwardDiff

# Constants
M = 1.0
a = 0.9 * M  # Spin parameter for Kerr black hole

# Initial conditions
r0 = 10.0       # Initial radial distance
θ0 = π / 2     # Equatorial plane
ϕ0 = 0.0       # Initial azimuthal angle
t0 = 0.0       # Initial time

# Metric function with corrected zeros
function metric(r, θ)
    Σ = r^2 + a^2 * cos(θ)^2
    Δ = r^2 - 2 * M * r + a^2
    sinθ = sin(θ)
    sin2θ = sinθ^2

    # Metric components
    g_tt = -(1 - (2 * M * r) / Σ)
    g_tϕ = - (2 * M * a * r * sin2θ) / Σ
    g_rr = Σ / Δ
    g_θθ = Σ
    g_ϕϕ = ( (r^2 + a^2)^2 - Δ * a^2 * sin2θ ) * sin2θ / Σ

    zero_value = zero(r)

    # Construct the metric tensor
    g = [
        g_tt        zero_value  zero_value  g_tϕ;
        zero_value  g_rr        zero_value  zero_value;
        zero_value  zero_value  g_θθ        zero_value;
        g_tϕ        zero_value  zero_value  g_ϕϕ
    ]
    return g
end

# Function to compute the Christoffel symbols using automatic differentiation
function compute_christoffel(r, θ)
    x = [r, θ]

    # Function to compute the metric components as a vector
    function g_vector(x)
        r, θ = x
        g = metric(r, θ)
        g_vec = []
        for α in 1:4
            for β in α:4  # Only upper triangle due to symmetry
                push!(g_vec, g[α, β])
            end
        end
        return g_vec
    end

    # Compute the Jacobian of the metric components with respect to coordinates
    g_vec = g_vector(x)
    dg_dx = ForwardDiff.jacobian(g_vector, x)

    # Reconstruct the metric tensor and its derivatives
    T = typeof(g_vec[1])
    g = Array{T}(undef, 4, 4)
    ∂g = zeros(T, 4, 4, 4)  # Partial derivatives with respect to t, r, θ, ϕ

    idx = 1
    for α in 1:4
        for β in α:4
            g[α, β] = g_vec[idx]
            g[β, α] = g[α, β]  # Symmetry
            ∂g[α, β, 2] = dg_dx[idx, 1]  # ∂g_{αβ}/∂r
            ∂g[β, α, 2] = ∂g[α, β, 2]    # Symmetry
            ∂g[α, β, 3] = dg_dx[idx, 2]  # ∂g_{αβ}/∂θ
            ∂g[β, α, 3] = ∂g[α, β, 3]    # Symmetry
            idx += 1
        end
    end

    g_inv = inv(g)

    # Compute the Christoffel symbols Γ^μ_{νλ}
    Γ = zeros(T, 4, 4, 4)
    for μ in 1:4
        for ν in 1:4
            for λ in 1:4
                sum = zero(T)
                for σ in 1:4  # Sum over all coordinates
                    term1 = ∂g[λ, σ, ν]   # ∂g_{λσ}/∂x^ν
                    term2 = ∂g[ν, σ, λ]   # ∂g_{νσ}/∂x^λ
                    term3 = ∂g[ν, λ, σ]   # ∂g_{νλ}/∂x^σ
                    sum += g_inv[μ, σ] * (term1 + term2 - term3)
                end
                Γ[μ, ν, λ] = 0.5 * sum
            end
        end
    end

    return Γ
end

# Set initial velocity components
v_r = -0.1     # Small inward radial velocity
v_θ = 0.0      # Polar velocity
v_ϕ = 0.5      # Azimuthal (angular) velocity

# Compute the metric at the initial position
g0 = metric(r0, θ0)

# Extract metric components
g_tt = g0[1, 1]
g_tϕ = g0[1, 4]
g_rr = g0[2, 2]
g_θθ = g0[3, 3]
g_ϕϕ = g0[4, 4]

# Compute coefficients for quadratic equation in v_t
A = g_tt
B = 2 * g_tϕ * v_ϕ
C = g_rr * v_r^2 + g_θθ * v_θ^2 + g_ϕϕ * v_ϕ^2

# Solve quadratic equation A * v_t^2 + B * v_t + C = 0
Δ = B^2 - 4 * A * C
if Δ < 0
    error("No real solution for v_t. Check initial conditions.")
end

v_t = (-B - sqrt(Δ)) / (2 * A)  # Choose the negative root for future-directed motion

# Complete the initial velocity vector
v = [v_t; v_r; v_θ; v_ϕ]

# Verify the normalization condition (should be close to zero)
norm = g_tt * v[1]^2 + 2 * g_tϕ * v[1] * v[4] + g_rr * v[2]^2 +
       g_θθ * v[3]^2 + g_ϕϕ * v[4]^2
println("Null normalization check: ", norm)

# Initial position and velocity vector
u0 = [t0, r0, θ0, ϕ0, v[1], v[2], v[3], v[4]]

# Event horizon radius for the Kerr metric
r_plus = M + sqrt(M^2 - a^2)  # Outer event horizon

# Condition function for the callback
function condition(u, t, integrator)
    r = u[2]
    return r - r_plus
end

# Affect function to terminate integration
function affect!(integrator)
    terminate!(integrator)
end

# Create the ContinuousCallback
termination_cb = ContinuousCallback(condition, affect!)

# Define the ODE function with corrections
function intprob!(du, u, p, λ)
    x = u[1:4]  # Position vector
    v = u[5:8]  # Velocity vector

    r_val, θ_val = x[2], x[3]

    du .= zero(eltype(u))

    # Compute the Christoffel symbols at the current position
    Γ = compute_christoffel(r_val, θ_val)

    # Compute the derivatives of velocity
    dv = zeros(eltype(u), 4)
    for μ in 1:4
        sum_ = zero(eltype(u))
        for ν in 1:4
            for λ in 1:4
                sum_ += Γ[μ, ν, λ] * v[ν] * v[λ]
            end
        end
        dv[μ] = -sum_
    end

    # Assign derivatives to du
    du[1:4] = v
    du[5:8] = dv
end


# Define the ODE problem and solve
tspan = (0.0, 100.0)
prob = ODEProblem(intprob!, u0, tspan)
sol = solve(prob, Tsit5(), callback=termination_cb, abstol=1e-10, reltol=1e-10)

# Extract trajectory data
t_vals = sol.t
r_vals = [u[2] for u in sol.u]
θ_vals = [u[3] for u in sol.u]
ϕ_vals = [u[4] for u in sol.u]

# Plot the trajectory in polar coordinates
pl = plot(
    θ -> r_plus,
    0:0.01:2π,
    lw = 2,
    legend = true,
    proj = :polar,
    color = :black,
    ylim=(0.0, maximum(r_vals) + 1.0),
    label="Event Horizon"
)

plot!(ϕ_vals, r_vals, lw=2, label="Photon Path", color=:blue)
display(pl)
