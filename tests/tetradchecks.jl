# to go with pitchangle.jl

# for checks ################################
# contravariant p^μ components
p_t = - p_t_hit
p_r = - p_r_hit
p_θ = - p_θ_hit
p_ϕ = - p_ϕ_hit

# contravariant p^μ at origin point on the disc
kμ = [p_t, p_r, p_θ, p_ϕ]

# check for norm
H1 = g[1,1]*p_t^2 + g[2,2]*p_r^2 + g[3,3]*p_θ^2 + g[4,4]*p_ϕ^2 + 2*g[1,4]*p_t*p_ϕ

# covariant p_μ at origin point on the disc
p_low_t = g_tt * p_t + g_tϕ * p_ϕ
p_low_r = g_rr * p_r
p_low_θ = g_θθ * p_θ
p_low_ϕ = g_tϕ * p_t + g_ϕϕ * p_ϕ

H2 = gtt * p_low_t^2 + grr * p_low_r^2 + gθθ * p_low_θ^2 + gϕϕ * p_low_ϕ^2 + 2 * gtϕ * p_low_t * p_low_ϕ

######################################################

e = lnrf_tetrad(r_hit, θ_hit, a)

# Transform momentum to LNRF frame
p_lnrf = zeros(4)
for μ in 1:4
    for ν in 1:4
        p_lnrf[μ] += e[μ,ν] * kμ[ν]
    end
end

p_spatial = p_lnrf[2:4]

# p(r)
p_spatial[1]
# p(θ)
p_spatial[2]
# p(ϕ)
p_spatial[3]

# inverse tetrad
inv(e)

# check: does p_BL[μ] = kμ[μ] --- YES
p_BL = inv(e) * p_lnrf

p_BL[1]
kμ[1]
# r component
p_BL[2]
kμ[2]
# θ component
p_BL[3]
kμ[3]
p_BL[4]
kμ[4]

p_spatial[1] * inv(e)[2, 2]
p_spatial[2] * inv(e)[3,3]

p_spatial[1] * sqrt(g_rr)
p_spatial[2] * sqrt(g_θθ)

p_spatial[2]
kμ[3] * sqrt(gθθ)

Q = p_spatial[3]^2 * (cos(μ_rad))^2

kμ[3]^2 * gθθ * (sin(μ_rad))^2 - (cos(μ_rad))^2 * kμ[2]^2 * grr

p_low_θ^2 * gθθ^3 * (sin(μ_rad))^2 - (cos(μ_rad))^2 * p_low_r^2 * grr^3

grr * p_low_r^2 + gθθ * p_low_θ^2 - 2 * gtϕ * E * L + gtt * E^2 + gϕϕ * L^2

g_rr * p_r^2 + g_θθ * p_low_θ^2


p_low_r
W = 2 * gtϕ * E * L - gtt * E^2 - gϕϕ * L^2
sqrt((W - gθθ * p_low_θ^2) / grr)

pθ = sqrt((Q + (cos(μ_rad))^2 * W * grr^2) / (gθθ * (gθθ^2 * (sin(μ_rad))^2 + (cos(μ_rad))^2 * grr^2)))

sqrt((W - gθθ * pθ^2)/(grr))