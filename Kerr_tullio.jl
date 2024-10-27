using Symbolics, Tullio, LinearAlgebra, DifferentialEquations, StaticArrays

# Define the metric as in the original code
@variables M t r θ ϕ
Σ = r^2 + a^2 * cos(θ)^2
delta = r^2 - 2*M*r + a^2

g = zeros(Num, (4, 4))
g[1, 1] = -(1 - 2 * M * r / Σ)
g[2, 2] = delta / Σ
g[3, 3] = Σ
g[4, 4] = (r^2 + a^2)^2 - delta * a^2 * sin(θ)^2
g[1, 4] = -2 * a * r * sin(θ)^2 / Σ
g[4, 1] = g[1, 4]

# Compute the Christoffel symbols
δ = map(Differential, [t, r, θ, ϕ])
function compute_christoffel(g, δ)
    ginv = inv(g)
    @tullio Γ[i, j, k] := 1 / 2 * ginv[i, m] * (
        expand_derivatives(
            δ[k](g[m, j]) + δ[j](g[m, k]) - δ[m](g[j, k])
        )
    )
    Γ
end
Γ = compute_christoffel(g, δ)

# Build the geodesic equation function
@variables v[1:4]
@tullio geodesic_eq[i] := -Γ[i, j, k] * v[j] * v[k]
geodesic_func = build_function(geodesic_eq, [t, r, θ, ϕ], v, M)

# Define ODE function using the built symbolic function
geo_func = eval(geodesic_func[1]) # Use out-of-place variant

function intprob!(du, u, params, λ)
    x = @view u[1:4]
    v = @view u[5:8]
    du[1:4] .= v
    du[5:8] .= geo_func(x, v, params...)
end

# Initial conditions
u0 = SVector(0.0, 200.0, π/2, 0.0, 
    normalize([0.0, -1.0, 0.0, 1.5e-4]))
λspan = (0.0, 500.0)

# Solve the ODE
prob = ODEProblem(intprob!, u0, λspan, [1.0])
sol = solve(prob, Tsit5())


