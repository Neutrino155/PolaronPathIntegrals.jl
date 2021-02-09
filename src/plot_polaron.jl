# plot_polaron.jl

# Interactive Figures
plotly()
Plots.PlotlyBackend()

# Static Figures
# pgfplots()

function plot_polaron(
    polaron,
    plot_type,
    xtickstep = 1,
    ytickstep = 1,
    ymax = 400.0,
    ylog = false,
)

    α = polaron.α
    v_values = polaron.v
    w_values = polaron.w
    Ω_range = polaron.Ω
    β_range = polaron.β
    T_range = polaron.T
    μ_values = polaron.μ
    Γ_values = polaron.Γ

    if plot_type == "mobility_frequency"
        x_values = Ω_range
        y_values = μ_values
        legend_values = T_range
        x_label = "Ω"
        y_label = "μ(Ω)"
        legend_label = "T = "
    elseif plot_type == "mobility_temperature"
        x_values = T_range
        y_values = reverse([[x[i] for x in μ_values] for i in 1:length(Ω_range)])
        legend_values = reverse(Ω_range)
        x_label = "T"
        y_label = "μ(T)"
        legend_label = "Ω = "
    elseif plot_type == "absorption_frequency"
        x_values = Ω_range
        y_values = Γ_values
        legend_values = T_range
        x_label = "Ω"
        y_label = "Γ(Ω)"
        legend_label = "T = "
    elseif plot_type == "absorption_temperature"
        x_values = T_range
        y_values = [[x[i] for x in Γ_values] for i in 1:length(Ω_range)]
        legend_values = Ω_range
        x_label = "T"
        y_label = "Γ(T)"
        legend_label = "Ω = "
    end

    if ylog
        plot = Plots.plot(
            x_values,
            y_values[1],
            xticks = round(x_values[1], digits = 1):xtickstep:round(x_values[end], digits = 1),
            yaxis = :log,
            xlabel = x_label,
            ylabel = y_label,
            label = legend_label * "$(round(legend_values[1], digits = 3))",
            legend = :bottomright,
            size = (1800, 1200),
            linewidth = 2,
            xtickfontsize = 18,
            ytickfontsize = 18,
            xguidefontsize = 20,
            yguidefontsize = 20,
            legendfontsize = 18,
            thickness_scaling = 1.2,
        )
        annotate!(plot,  Plots.xlims(plot)[1] + 0.8 * (Plots.xlims(plot)[2] - Plots.xlims(plot)[1]), 0.5e-2, ("\$\\alpha\$ = $α", :black, :right, 20))
        for i = 2:length(legend_values)
            Plots.plot!(
                plot,
                x_values,
                y_values[i],
                label = legend_label * "$(round(legend_values[i], digits = 3))",
                linewidth = 2,
            )
        end

    else
        plot = Plots.plot(
            x_values,
            y_values[1],
            xticks = round(x_values[1], digits = 1):xtickstep:round(x_values[end], digits = 1),
            ylims = (0.0, ymax),
            xlabel = x_label,
            ylabel = y_label,
            label = legend_label * "$(round(legend_values[1], digits = 3))",
            legend = :bottomright,
            size = (1800, 1200),
            linewidth = 2,
            xtickfontsize = 18,
            ytickfontsize = 18,
            xguidefontsize = 20,
            yguidefontsize = 20,
            legendfontsize = 18,
            thickness_scaling = 1.2,
        )
        annotate!(plot, Plots.xlims(plot)[1] + 0.4 * (Plots.xlims(plot)[2] - Plots.xlims(plot)[1]), Plots.ylims(plot)[1] + 0.95 * (Plots.ylims(plot)[2] - Plots.ylims(plot)[1]), ("\$\\alpha\$ = $(round(α, digits = 3))", :black, :right, 20))
        for i = 2:Int(length(legend_values))
            Plots.plot!(
                plot,
                x_values,
                y_values[i],
                label = legend_label * "$(round(legend_values[i], digits = 3))",
                linewidth = 2,
            )
        end
    end
    return plot
end
