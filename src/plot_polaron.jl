# plot_polaron.jl

# Interactive Figures
# plotly()
# Plots.PlotlyBackend()

# Static Figures
# pyplot()
# default(fmt=:svg)

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

function plot_polaron_interactive(ϵ_optic, ϵ_static, phonon_freq, m_eff; t = 100:15:400, f = 0:1:20)

    fig = GLMakie.Figure(resolution = (1920, 1080), backgroundcolor = RGBf0(0.7, 0.7, 0.7))

    energy_axis = fig[1, 3] = Axis(fig, xlabel = "T [K]", ylabel = "F(T) [meV]", title = "Free Energy")
    vary_axis = fig[2, 3] = Axis(fig, xlabel = "T [K]", ylabel = "v/w(T) [Hz]", title = "Variational Parameters")
    mass_axis = fig[1, 4] = Axis(fig, xlabel = "T [K]", ylabel = "M(T) [m_e kg]", title = "Fictitious Mass")
    spring_axis = fig[2, 4] = Axis(fig, xlabel = "T [K]", ylabel = "κ(T) [me kg/s^2]", title = "Fictitious Spring Constant")
    mobility_axis_temp = fig[1, 2] = Axis(fig, xlabel = "T [K]", ylabel = "μ(T) [cm^2/Vs]", title = "Mobility vs Temperature")
    absorption_axis_temp = fig[2, 2] = Axis(fig, xlabel = "T [K]", ylabel = "Γ(T) [cm^-1]", title = "Optical Absorption vs Temperature")
    mobility_axis_freq = fig[1, 5] = Axis(fig, xlabel = "Ω [Hz]", ylabel = "μ(Ω) [cm^2/Vs]", title = "Mobility vs Frequency")
    absorption_axis_freq = fig[2, 5] = Axis(fig, xlabel = "Ω [Hz]", ylabel = "Γ(Ω) [cm^-1]", title = "Optical Absorption vs Frequency")

    lsgrid = labelslidergrid!(
        fig,
        ["ϵ_0", "ϵ_∞", "m", "ω"],
        [1:0.01:100, 1:0.01:100, 0.1:0.01:1, 0.1:0.01:10],
        formats = [x -> "$(round(x, digits = 3)) $s" for s in ["F/m", "F/m", "m_e", "THz"]],
        width = 1000
    )

    set_close_to!(lsgrid.sliders[1], ϵ_static)
    set_close_to!(lsgrid.sliders[2], ϵ_optic)
    set_close_to!(lsgrid.sliders[3], m_eff)
    set_close_to!(lsgrid.sliders[4], phonon_freq)

    fig[3, 1:5] = lsgrid.layout

    Ω_slider = Slider(fig[1:2, 1], range = f, startvalue = 0.0, horizontal = false, tellwidth = true, height = nothing, width = Auto())
    freq_sliders = [s.value for s in vcat(lsgrid.sliders , [Ω_slider])]
    Ω_slider_label = Label(fig[1:2, 0], lift(s -> "Ω: $(s[]) Hz", Ω_slider.value), rotation = pi/2)

    polaron_temp = lift(freq_sliders...) do slvalues...
        make_polaron(slvalues[2], slvalues[1], slvalues[4] * 1e12, slvalues[3]; temp = t, efield_freq = slvalues[5], verbose = false)
    end

    T = lift(p -> p.T, polaron_temp)
    v = lift(p -> Vector{Float64}(p.v), polaron_temp)
    w = lift(p -> Vector{Float64}(p.w), polaron_temp)
    κ = lift(p -> Vector{Float64}(p.κ), polaron_temp)
    M = lift(p -> Vector{Float64}(p.M), polaron_temp)
    F = lift(p -> Vector{Float64}(p.F), polaron_temp)
    μ_T = lift(p -> vcat.(p.μ...)[1], polaron_temp)
    Γ_T = lift(p -> vcat.(p.Γ...)[1], polaron_temp)

    T_slider = Slider(fig[1:2, 7], range = t, startvalue = 300, horizontal = false, tellwidth = true, height = nothing, width = Auto())
    T_slider_label = Label(fig[1:2, 8], lift(s -> "T: $(s[]) K", T_slider.value), rotation = -pi/2)

    temp_sliders = [s.value for s in vcat(lsgrid.sliders , [T_slider])]

    polaron_freq = lift(temp_sliders...) do slvalues...
        make_polaron(slvalues[2], slvalues[1], slvalues[4] * 1e12, slvalues[3]; temp = slvalues[5], efield_freq = f, verbose = false)
    end

    Ω = lift(p -> p.Ω, polaron_freq)
    μ_Ω = lift(p -> Vector{Float64}(p.μ[1]), polaron_freq)
    Γ_Ω = lift(p -> Vector{Float64}(p.Γ[1]), polaron_freq)

    scatterlines!(mobility_axis_temp, T, μ_T)
    scatterlines!(absorption_axis_temp, T, Γ_T)
    scatterlines!(energy_axis, T, F)
    v_plot = scatterlines!(vary_axis, T, v, label = "v")
    w_plot = scatterlines!(vary_axis, T, w, label = "w", marker = :diamond)
    scatterlines!(spring_axis, T, κ)
    scatterlines!(mass_axis, T, M)

    axislegend(vary_axis, [v_plot, w_plot], ["v", "w"], position = :rb)

    on(polaron_temp) do x
        autolimits!(mobility_axis_temp)
        autolimits!(absorption_axis_temp)
        autolimits!(energy_axis)
        autolimits!(vary_axis)
        autolimits!(spring_axis)
        autolimits!(mass_axis)
    end

    scatterlines!(mobility_axis_freq, Ω, μ_Ω)
    scatterlines!(absorption_axis_freq, Ω, Γ_Ω)

    on(polaron_freq) do x
        autolimits!(mobility_axis_freq)
        autolimits!(absorption_axis_freq)
    end

    fig
end
