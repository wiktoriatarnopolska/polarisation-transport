function disc_hits_plot(r_in, r_out, r_horizon, x_hits, y_hits, x_vals, y_vals)

    θ_values = range(0, 2π, length=500)
    x_inner = [r_in * cos(θ) for θ in θ_values]
    y_inner = [r_in * sin(θ) for θ in θ_values]
    x_outer = [r_out * cos(θ) for θ in θ_values]
    y_outer = [r_out * sin(θ) for θ in θ_values]

    pl_disc = plot(
    x_outer, y_outer,
    seriestype = :shape,
    fillcolor = :mediumorchid2,
    linecolor = :rebeccapurple,
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    title = "Disc Hits Visualization",
    label = "Disc",
    xlim = (-20, 20),
    ylim = (-20, 20),
    )

    # Fill the inner disc area
    plot!(
    x_inner, y_inner,
    seriestype = :shape,
    fillcolor = :white,
    linecolor = :rebeccapurple,
    label = false
    )

    # Overlay the event horizon
    circle_x = [r_horizon * cos(θ) for θ in 0:0.01:2π]
    circle_y = [r_horizon * sin(θ) for θ in 0:0.01:2π]
    plot!(
    circle_x, circle_y,
    lw = 2,
    color = :purple4,
    label = "Event horizon"
    )

    # Overlay the ISCO
    isco_x = [r_in * cos(θ) for θ in 0:0.01:2π]
    isco_y = [r_in * sin(θ) for θ in 0:0.01:2π]
    plot!(
    isco_x, isco_y,
    lw = 2,
    color = :mediumslateblue,
    label = "ISCO"
    )

    # Overlay the hit points
    scatter!(
    pl_disc,
    x_hits, y_hits,
    color = :tan1,
    markerstrokecolor = :black,
    markersize = 6,
    label = "Disc Hits",
    )

    # Plot the trajectory
    plot!(
    pl_disc,
    x_vals, y_vals,
    lw = 2.0,
    linecolor = :steelblue,
    label = "incident ray"
    )

end

function conservations_plot(sol, E_vals, L_vals, Q_vals, H_vals)
# Extract affine parameter values
    λ_vals = sol.t

    abs_E_diff = abs.(E_vals .- E_vals[1])
    abs_E_diff[abs_E_diff .== 0] .= 1e-10  # Replace zeros with small values

    abs_L_diff = abs.(L_vals .- L_vals[1])
    abs_L_diff[abs_L_diff .== 0] .= 1e-10  # Replace zeros with small values

    abs_Q_diff = abs.(Q_vals .- Q_vals[1])
    abs_Q_diff[abs_Q_diff .== 0] .= 1e-10  # Replace zeros with small values

    abs_H_diff = abs.(H_vals .- H_vals[1])
    abs_H_diff[abs_H_diff .== 0] .= 1e-10  # Replace zeros with small values

    pl_conserved = plot(
        xlabel="Affine Parameter λ",
        ylabel="|ΔQuantity| (log scale)",
        title="Quantity Difference",
        legend=:outerright,
        yscale=:log10,
    )

    plot!(
        pl_conserved,
        λ_vals, abs_E_diff,
        label="|ΔE|",
        lw=1.5
    )

    # pl_Lconserved = plot(
    #     xlabel="Affine Parameter λ",
    #     ylabel="|ΔL| (log scale)",
    #     title="Momentum Difference |L₀ - L(λ)|",
    #     legend=:outerright,
    #     yscale=:log10,
    # )

    plot!(
        pl_conserved,
        λ_vals, abs_L_diff,
        label="|ΔL|",
        lw=1.5
    )

    # pl_Qconserved = plot(
    #     xlabel="Affine Parameter λ",
    #     ylabel="|ΔQ| (log scale)",
    #     title="Carter Const. Difference |Q₀ - Q(λ)|",
    #     legend=:outerright,
    #     yscale=:log10,
    # )

    plot!(
        pl_conserved,
        λ_vals, abs_Q_diff,
        label="|ΔQ|",
        lw=1.5
    )

    # pl_Hconserved = plot(
    #     xlabel = "Affine Parameter λ",
    #     ylabel = "|ΔH| (log scale)",
    #     title = "Hamiltonian Difference |H₀ - H(λ)|",
    #     legend = :outerright,
    #     yscale = :log10
    # )

    plot!(
            pl_conserved,
            λ_vals, abs_H_diff,
            label = "|ΔH|",
            lw = 1.5
        )
end


export disc_hits_plot, conservations_plot