# plot_polaron.jl

plotly()
Plots.PlotlyBackend()

function plot_polaron(polaron)
    Ω = polaron.Ω
    T = polaron.T
    if typeof(Ω) != Float64
        mobility_plot = Plots.plot(Ω, polaron.μ, title = "Polaron Mobility", xlabel = "Ω", ylabel = "μ(Ω)")
        absorption_plot = Plots.plot(Ω, polaron.Γ, title = "Polaron Optical Absorption", xlabel = "Ω", ylabel = "Γ(Ω)")
        Plots.plot(mobility_plot, absorption_plot, layout = (1, 2), legend = false)
    elseif typeof(T) != Float64
        mobility_plot = Plots.plot(T, polaron.μ, title = "Polaron Mobility", xlabel = "T", ylabel = "μ(T)")
        absorption_plot = Plots.plot(T, polaron.Γ, title = "Polaron Optical Absorption", xlabel = "T", ylabel = "Γ(T)")
        Plots.plot(mobility_plot, absorption_plot, layout = (1, 2), legend = false)
    end

end
