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
    Z = polaron.Z # Mobilities
    σ = polaron.σ # Optical absorptions

    # Plot Mobility versus Frequency for N different Temperatures. 
    σ_Ω = Plots.plot(Ω, [real(σ[:, i]) for i in 1:Int(floor(length(T)/(N))):length(T)], label = hcat(["T = $i" for i in T][1:Int(floor(length(T)/(N))):end]...), title = "Conductivity σ(Ω)", xlabel = "Ω / ω", ylabel = "σ(Ω)", legend = true, minorgrid = true, linewidth = 1.8, xtickfontsize = 15, ytickfontsize = 15, xguidefontsize = 15, yguidefontsize = 15, legendfontsize = 15, thickness_scaling = 1.5, size = (600, 600), linestyle  = :solid)
    cur_colors = theme_palette(:default)
    for i in 1:Int(floor(length(T)/(N))):length(T)
        Plots.plot!(σ_Ω, Ω, imag(σ[:, i]), label = false, linestyle = :dash, color = cur_colors[Int(i-1+Int(floor(length(T)/(N))))÷Int(floor(length(T)/(N)))])
    end

    # Plot Mobility versus Temperature for N different Frequencies.
    σ_T = Plots.plot(T, [real.(σ[i, :]) for i in 1:Int(floor(length(Ω)/(N))):length(Ω)], label = hcat(["Ω = $i" for i in Ω][1:Int(floor(length(Ω)/(N))):end]...), title = "Conductivity σ(T)",  xlabel = "T / ω", ylabel = "σ(T)", legend = true, minorgrid = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))
    cur_colors = theme_palette(:default)
    for i in 1:Int(floor(length(Ω)/(N))):length(Ω)
        Plots.plot!(σ_T, T, imag(σ[i, :]), label = false, linestyle = :dash, color = cur_colors[Int(i-1+Int(floor(length(Ω)/(N))))÷Int(floor(length(Ω)/(N)))])
    end

    # Plot Asbsorption versus Frequency for N different Temperatures.
    Z_Ω = Plots.plot(Ω, [real(Z[:, i]) for i in 1:Int(floor(length(T)/(N))):length(T)], label = hcat(["T = $i" for i in T][1:Int(floor(length(T)/(N))):end]...), title = "Impedence Z(Ω)",  xlabel = "Ω / ω", ylabel = "Z(Ω)", legend = true, minorgrid = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))
    cur_colors = theme_palette(:default)
    for i in 1:Int(floor(length(T)/(N))):length(T)
        Plots.plot!(Z_Ω, Ω, imag(Z[:, i]), label = false, linestyle = :dash, color = cur_colors[Int(i-1+Int(floor(length(T)/(N))))÷Int(floor(length(T)/(N)))])
    end
    # Plot Absorption versus Temperature for N different Frequencies.
    Z_T = Plots.plot(T, [real(Z[i, :]) for i in 1:Int(floor(length(Ω)/(N))):length(Ω)], label = hcat(["Ω = $i" for i in Ω][1:Int(floor(length(Ω)/(N))):end]...), title = "Impedence Z(T)",  xlabel = "T / ω", ylabel = "Z(T)", legend = true, minorgrid = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))
    cur_colors = theme_palette(:default)
    for i in 1:Int(floor(length(Ω)/(N))):length(Ω)
        Plots.plot!(Z_T, T, imag(Z[i, :]), label = false, linestyle = :dash, color = cur_colors[Int(i-1+Int(floor(length(Ω)/(N))))÷Int(floor(length(Ω)/(N)))])
    end

    # Plot Spring Constant versus Temperature.
    κ_T = Plots.plot(T, κ, title = "Spring Constant κ(T)", xlabel = "T / ω", ylabel = "κ / m_e", legend = false, minorgrid = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Fictitious Mass versus Temperature.
    M_T = Plots.plot(T, M, title = "Fictitious Mass M(T)", xlabel = "T / ω", ylabel = "M(T)", legend = false, minorgrid = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Variational Parameters v & w versus Temperature.
    vw_T = Plots.plot(T, [v, w], label = ["v" "w"], title = "Variational Parameters v(T) & w(T)", xlabel = "T(K)", ylabel = "v & w (s^{-1})", legend = true, minorgrid = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Plot Free Energy versus Temperature.
    F_T = Plots.plot(T, F, title = "Free Energy F(T)", xlabel = "T(K)", ylabel = "F(T) (meV)", legend = false, minorgrid = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, thickness_scaling = 1.2, size = (600, 600))

    # Combine all plots above as subplots for combined view.
    all_plots = Plots.plot(σ_Ω, σ_T, Z_Ω, Z_T, κ_T, M_T, vw_T, F_T, layout = (2, 4), size = (1800, 900), minorgrid = true, linewidth = 1.5, xtickfontsize = 10, ytickfontsize = 8, xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 10, thickness_scaling = 1.2)

    # Return all plots.
    all_plots, σ_Ω, σ_T, Z_Ω, Z_T, κ_T, M_T, vw_T, F_T
end

"""
save_polaron_plots(plots::Array{Plot}, path::String, ext::String)

    Saves all the plots created from plot_polaron() at the specified path with the specified extension.
"""

function save_polaron_plots(plots, path, ext = "svg")
    for plot in plots
        Plots.savefig(plot, path * "$(plot)." * ext)
    end
end

function plot_polaron_interactive(ϵ_optic, ϵ_static, phonon_freq, m_eff; t = 100:15:400, f = 0:1:20)

    fig = GLMakie.Figure(resolution = (1920, 1080), backgroundcolor = RGBf0(0.9, 0.9, 0.9))

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
        [1:0.01:100, 1:0.01:100, 0.1:0.01:1, 0.01:0.01:100],
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
