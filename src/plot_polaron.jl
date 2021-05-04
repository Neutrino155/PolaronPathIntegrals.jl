# plot_polaron.jl

# Interactive Figures
# plotly()
# Plots.PlotlyBackend()

# Static Figures
pyplot()
default(fmt=:svg)

"""
plot_polaron(polaron::Polaron; N::Int)

    Creates the following plots of a given Polaron mutable struct:

        • μ_Ω: Mobility μ(Ω) (cm^2/Vs) versus Frequency Ω/ω (multiple of phonon frequency). Includes plots at N different evenly-spaced Temperatures T (K).
        • μ_T: Mobility μ(T) (cm^2/Vs) versus Temperature T (K). Includes plots at N different evenly-spaced Frequencies Ω/ω (multiple of phonon frequency).
        • Γ_Ω: Optical Absorption Γ(Ω) (cm^-1) versus Frequency Ω/ω (multiple of phonon frequency). Includes plots at N different evenly-spaced Temperatures T (K).
        • Γ_T: Optical Absorption Γ(Ω) (cm^-1) versus Temperature T (K). Includes plots at N different evenly-spaced Frequencies Ω/ω (multiple of phonon frequency).
        • κ_T: Spring Constant κ(T)/m_e (kg/s^2) (multiple of electron masses) versus Temperature T (K).
        • M_T: Fictitious mass M(T)/m_e (kg) (multiple of electron masses) versus Temperature T (K).
        • vw_T: Variational Parameters v and w (s^-1) versus Temperature T (K).
        • F_T: Free Energy (meV) versus Temperature T (K).
        ⦿ all_plots: A final plot that includes all the above as different subplots.

    polaron is a Polaron Mutable Stuct type created from make_polaron(). N is an Integer that specifies how many different temperatures or frequencies to plots on the Mobility and Optical Absorption plots.

    returns plots: all_plots, μ_Ω, μ_T, Γ_Ω, Γ_T, κ_T, M_T, vw_T, F_T
"""

function plot_polaron(polaron; N = 6)

    # Extract polaron data.
    α = polaron.α # Alpha parameter
    v = polaron.v # v variational parameter
    w = polaron.w # w variational parameter
    F = polaron.F # Free energy
    κ = polaron.κ # Spring constant
    M = polaron.M # Fictitious mass
    Ω = polaron.Ω # Electric field frequencies
    β = polaron.β # Reduced thermodynamic betas
    T = polaron.T # Temperatures
    μ = polaron.μ # Mobilities
    Γ = polaron.Γ # Optical absorptions

    # Plot Mobility versus Frequency for N different Temperatures.
    μ_Ω = Plots.plot(Ω, μ[1:Int(floor(length(T)/(N - 1))):end], label = hcat(["T = $i" for i in T][1:Int(floor(length(T)/(N - 1))):end]...), title = "Mobility \$\\mu(\\Omega)\$", ylim = (0, maximum(μ)[end] * 1.1), xlabel = "\$\\Omega\$ / \$\\omega\$", ylabel = "\$\\mu(\\Omega)\$ \$(cm^2/Vs)\$", legend = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Mobility versus Temperature for N different Frequencies.
    μ_T = Plots.plot(T, vcat.(μ...)[1:Int(floor(length(Ω)/(N - 1))):end], label = hcat(["Ω = $i" for i in Ω][1:Int(floor(length(Ω)/(N - 1))):end]...), title = "Mobility \$\\mu(T)\$", ylim = (0, maximum(μ)[end] * 1.1),  xlabel = "\$T\$ \$(K)\$", ylabel = "\$\\mu(T)\$ \$(cm^2/Vs)\$", legend = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Asbsorption versus Frequency for N different Temperatures.
    Γ_Ω = Plots.plot(Ω, Γ[1:Int(floor(length(T)/(N - 1))):end], label = hcat(["T = $i" for i in T][1:Int(floor(length(T)/(N - 1))):end]...), title = "Optical Absorption \$\\Gamma(\\Omega)\$", ylim = (0, sort(maximum(vcat.(Γ...)))[end - 1] * 1.1),  xlabel = "\$\\Omega\$ / \$\\omega\$", ylabel = "\$\\Gamma(\\Omega)\$ \$(cm^{-1})\$", legend = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Absorption versus Temperature for N different Frequencies.
    Γ_T = Plots.plot(T, vcat.(Γ...)[1:Int(floor(length(Ω)/(N - 1))):end], label = hcat(["Ω = $i" for i in Ω][1:Int(floor(length(Ω)/(N - 1))):end]...), title = "Optical Absorption \$\\Gamma(T)\$", ylim = (0, sort(maximum(vcat.(Γ...)))[end - 1] * 1.1),  xlabel = "\$T\$ \$(K)\$", ylabel = "\$\\Gamma(T)\$ \$(cm^{-1})\$", legend = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Spring Constant versus Temperature.
    κ_T = Plots.plot(T, κ, title = "Spring Constant \$\\kappa(T)\$", xlabel = "\$T\$ \$(K)\$", ylabel = "\$\\kappa(T)\$ / \$m_e\$ \$(kg/s^2)\$", legend = false, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Fictitious Mass versus Temperature.
    M_T = Plots.plot(T, M, title = "Fictitious Mass \$M(T)\$", xlabel = "\$T\$ \$(K)\$", ylabel = "\$M(T)\$ / \$m_e\$ \$(kg)\$", legend = false, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Variational Parameters v & w versus Temperature.
    vw_T = Plots.plot(T, [v, w], label = ["v" "w"], title = "Variational Parameters \$v(T)\$ & \$w(T)\$", xlabel = "\$T\$ \$(K)\$", ylabel = "\$v\$ & \$w\$ \$(s^{-1})\$", legend = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Free Energy versus Temperature.
    F_T = Plots.plot(T, F, title = "Free Energy \$F(T)\$", xlabel = "\$T\$ \$(K)\$", ylabel = "\$F(T)\$ \$(meV)\$", legend = false, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Combine all plots above as subplots for combined view.
    all_plots = Plots.plot(μ_Ω, μ_T, Γ_Ω, Γ_T, κ_T, M_T, vw_T, F_T, layout = (2, 4), size = (1800, 900), linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 8, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 10, thickness_scaling = 1.2)

    # Return all plots.
    all_plots, μ_Ω, μ_T, Γ_Ω, Γ_T, κ_T, M_T, vw_T, F_T
end

"""
save_polaron_plots(plots::Array{Plot}, path::String, ext::String)

    Saves all the plots created from plot_polaron() at the specified path with the specified extension.
"""

function save_polaron_plots(plots, path, ext = "svg")

    for plot in plots
        savefig(plot, path * "$(plot)." * ext)
    end
end

# function plot_polaron_interactive(ϵ_optic, ϵ_static, phonon_freq, m_eff)
#
#     p = (:ϵ_∞, :ϵ_0, :ω, :m)
#     i = (4.5, 24.1, 2.25, 0.12)
#     pl = ((0.1, 10.0), (0.1, 10.0), (0.1, 10.0), (0.1, 1.0))
#     al = ((0, 300), (0, 10), (0, 1000))
#     a = Dict(
#         "ℜ_title" => "Mobility",
#         "ℜ_x" => "μ(T)",
#         "ℜ_y" => "μ(Ω)",
#         "ℜ_z" => "μ(T, Ω)",
#         "ℑ_title" => "Absorption",
#         "ℑ_x" => "Γ(T)",
#         "ℑ_y" => "Γ(Ω)",
#         "ℑ_z" => "Γ(T, Ω)",
#     )
#
#     polaron(T, Ω, ϵ_optic, ϵ_static, phonon_freq, m_eff) = make_polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff; temp = T, efield_freq = Ω)
#
#     vf = viewfunction(polaron.μ, polaron.Γ, p, i, pl, a, axis_limits = al, len = 1000)
# end
