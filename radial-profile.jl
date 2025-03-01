using Plots

rout = 10.0
###### a = 0.998 KERR CASE ########################################################################

# Calculate ISCO (in units of r_g)
am = 0.998 # Angular momentum a/M
# 2 digits (1.24)
isco = isco_radius(am)
println("ISCO radius: $isco r_g")

# Generate radii from ISCO to 1000 r_g (radii are now dimensionless, in units of r_g)
radii = range(isco, stop=rout, length=100)  # Radii in units of r_g

# Calculate Novikov-Thorne radial profile for each radius
profiles = [novikov_thorne_profile(r, isco, M, am, M_dot, f_col) for r in radii]

# Plotting the radial profile
plot(radii, profiles, xlabel=L"\textrm{Radius} \, (\textrm{r_g})", ylabel=L"\textrm{Temperature \, (keV)}", title=L"\textrm{Novikov-Thorne \, \,   Radial \, \,  Profile}",
label = "a = $am", color=:darkred, linewidth = 2)
# Marking the ISCO on the plot
vline!([isco], linestyle=:dash, color=:darkred, label = "ISCO")
#savefig("radial_profile.png")

###### a = 0.9 KERR CASE ###########################

am2 = 0.9
# Calculate ISCO (in units of r_g)
isco2 = isco_radius(am2)
println("ISCO radius: $isco2 r_g")
# Generate radii from ISCO to 1000 r_g (radii are now dimensionless, in units of r_g)
radii2 = range(isco2, stop=rout, length=100)  # Radii in units of r_g
# Calculate Novikov-Thorne radial profile for each radius
profiles2 = [novikov_thorne_profile(r, isco2, M, am2, M_dot, f_col) for r in radii2]
# Plotting the radial profile
plot!(radii2, profiles2, label = "a = $am2", color=:goldenrod1, linewidth=2)
# Marking the ISCO on the plot
vline!([isco2], linestyle=:dash, color=:goldenrod1, label = "ISCO")
# savefig("radial_profile.png")

###### a = 0.0 SCHWARZSCHILD CASE ###########################

am1 = 0.0
# Calculate ISCO (in units of r_g)
isco1 = isco_radius(am1)
println("ISCO radius: $isco1 r_g")
# Generate radii from ISCO to 1000 r_g (radii are now dimensionless, in units of r_g)
radii1 = range(isco1, stop=rout, length=100)  # Radii in units of r_g
# Calculate Novikov-Thorne radial profile for each radius
profiles1 = [novikov_thorne_profile(r, isco1, M, am1, M_dot, f_col) for r in radii1]
# Plotting the radial profile
plot!(radii1, profiles1, label = "a = $am1", color=:royalblue1, linewidth=2, frame=:box)
# Marking the ISCO on the plot
vline!([isco1], linestyle=:dash, color=:royalblue1, label = "ISCO")
# savefig("radial_profile.png")



# ###### LOG SCALE ###################################################
# # Plotting the radial profile with logarithmic scaling
# plot(radii, profiles, xlabel="Radius (r_g)", ylabel="Temperature (keV)", 
#      title="Novikov-Thorne Radial Profile", xscale=:log10, yscale=:log10, label = "a = 0.998", colour=:red)
# plot!(radii1, profiles1, label = "a = 0.0", color=:blue)
# plot!(radii2, profiles2, label = "a = 0.9", color=:green)
# # Marking the ISCO on the plot
#  vline!([isco], label="ISCO", linestyle=:dash, color=:red)
#  vline!([isco1], label="ISCO", linestyle=:dash, color=:blue)
#  vline!([isco2], label="ISCO", linestyle=:dash, color=:green)
# # savefig("radial_profile_logscale.png")

