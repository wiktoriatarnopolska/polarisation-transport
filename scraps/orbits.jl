using Plots
using LaTeXStrings

spins = collect(-1.0:0.001:1.0)

event_horizons = []
iscos = []
photon_circ_orbit = []
marginally_bound_orb = []

for s in spins
    r_eh = horizon(s)
    r_isco = isco_radius(s)
    r_ph = 2 * (1 + cos((2/3) * acos(-s)))
    r_mb = 2 - s + 2 * sqrt(1 - s)

    push!(event_horizons, r_eh)
    push!(iscos, r_isco)
    push!(photon_circ_orbit, r_ph)
    push!(marginally_bound_orb, r_mb)
end

p = plot(
    spins, event_horizons,
    label = L"\textrm{Event \, horizon}",
    xlabel = L"\mathrm{Spin} \, a \, (\mathrm{dimensionless})",
    ylabel = L"\mathrm{Radius} \, (r_\textrm{g})",
    legend=:topright,
    linecolor=:black,
    linewidth = 2,
    frame=:box
)

plot!(
    p, spins, iscos,
    label = L"\textrm{ISCO}",
    linecolor=:darkred,
    linewidth = 2,
    linestyle=:dash
)

plot!(
    p, spins, photon_circ_orbit,
    label = L"\textrm{Unstable \, photon \, orbit}",
    linecolor=:blueviolet,
    linewidth = 2,
    linestyle =:dot
)

