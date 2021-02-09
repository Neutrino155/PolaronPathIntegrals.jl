# FHIP_mobility.jl

"""
Reproduces figures 1-3 in the 1962 FHIP paper, and extends them to a larger range of frequencies and temperatures.
"""

using PolaronPathIntegrals
using Plots
using SpecialFunctions
using LaTeXStrings
using DSP
using Images

# Interactive Figures
plotly()
Plots.PlotlyBackend()

# Static Figures
# pgfplots()

function plot_figure(
    Ω_range,
    α,
    w,
    v;
    β_range = [
        gamma(i) for i in range(1.461632144968, stop = 5.89251869634377, length = 7)
    ],
    xtickstep = 1,
    ytickstep = 1,
    ylog = false,
    zero = false,
)

    if ylog
        zero = false
    end

    if !zero
        # Imχs = [[PolaronPathIntegrals.ℑχ(Ω, β, α, v, w) for Ω in Ω_range] for β in β_range]
        Imχs = [imag(hilbert(vcat(repeat([Float64(PolaronPathIntegrals.ℑχ(Ω_range[1], β, α, v, w))], 2 * length(Ω_range) + 1), [Float64(PolaronPathIntegrals.ℑχ(Ω, β, α, v, w)) for Ω in Ω_range], repeat([Float64(PolaronPathIntegrals.ℑχ(Ω_range[end], β, α, v, w))], 2 * length(Ω_range) + 1))))[2 * length(Ω_range) + 2:end - 2 * length(Ω_range) - 1] for β in β_range]
        @show(Imχs[1])
    else
        Imχs = [vcat([0.0], [PolaronPathIntegrals.ℑχ(Ω, β, α, v, w) for Ω in Ω_range]) for β in β_range]
        Ω_range = vcat([0.0], Ω_range)
    end

    if ylog
        plot = Plots.plot(
            Ω_range,
            Imχs[1],
            xticks = Ω_range[1]:xtickstep:Ω_range[end],
            yaxis = :log,
            xlabel = "Ω",
            ylabel = "Imχ(Ω)",
            label = "\$\\beta\$ = $(round(β_range[1], digits = 3))",
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
        annotate!(plot,  Plots.xlims(plot)[1] + 0.8 * (Plots.xlims(plot)[2] - Plots.xlims(plot)[1]), 0.5e-2, ("\$\\alpha\$ = $α, v = $v, w = $w", :black, :right, 20))
        for i = 2:length(Imχs)
            Plots.plot!(
                plot,
                Ω_range,
                Imχs[i],
                label = "\$\\beta\$ = $(round(β_range[i], digits = 3))",
                linewidth = 2,
            )
        end

    else
        plot = Plots.plot(
            Ω_range,
            Imχs[1],
            xticks = Ω_range[1]:xtickstep:Ω_range[end],
            yticks = floor(min((Imχs...)...)):ytickstep:ceil(max((Imχs...)...)),
            xlabel = "Ω",
            ylabel = "Imχ(Ω)",
            label = "\$\\beta\$ = $(round(β_range[1], digits = 3))",
            legend = :topleft,
            size = (1800, 1200),
            linewidth = 2,
            xtickfontsize = 18,
            ytickfontsize = 18,
            xguidefontsize = 20,
            yguidefontsize = 20,
            legendfontsize = 18,
            thickness_scaling = 1.2,
        )
        annotate!(plot, Plots.xlims(plot)[1] + 0.4 * (Plots.xlims(plot)[2] - Plots.xlims(plot)[1]), Plots.ylims(plot)[1] + 0.95 * (Plots.ylims(plot)[2] - Plots.ylims(plot)[1]), ("\$\\alpha\$ = $α, v = $v, w = $w", :black, :right, 20))
        for i = 2:length(Imχs)
            Plots.plot!(
                plot,
                Ω_range,
                Imχs[i],
                label = "\$\\beta\$ = $(round(β_range[i], digits = 3))",
                linewidth = 2,
            )
        end
    end
    return plot
end

"""
Replicated Plots
"""

setprecision(BigFloat, 64)

# replicated_figure_one = plot_figure(
#     range(3.0, stop = 9.0, length = 1000),
#     3.0,
#     2.5,
#     3.4;
#     β_range = [100.0],
# )
#
# replicated_figure_two = plot_figure(
#     range(2.0, stop = 22.0, length = 1000),
#     5.0,
#     2.1,
#     4.0;
#     β_range = [100.0],
#     xtickstep = 2,
#     ytickstep = 2,
# )
#
# replicated_figure_three = plot_figure(
#     range(6.0, stop = 28.0, length = 1000),
#     7.0,
#     1.6,
#     5.8;
#     β_range = [100.0],
#     xtickstep = 2,
#     ylog = true,
# )

"""
Extended Plots
"""

# extended_figure_one = plot_figure(
#     range(3.0, stop = 9.0, length = 1000),
#     3.0,
#     2.5,
#     3.4;
#     ylog = false,
#     xtickstep = 1,
#     ytickstep = 1,
#     zero = false,
# )

extended_figure_two = plot_figure(
    range(2.0, stop = 28.0, length = 1000),
    5.0,
    2.1,
    4.0;
    ylog = false,
    xtickstep = 2,
    ytickstep = 5,
    zero = false,
)
#
# extended_figure_three = plot_figure(
#     range(6.0, stop = 28.0, length = 1000),
#     7.0,
#     1.6,
#     5.8;
#     ylog = false,
#     xtickstep = 2,
#     ytickstep = 5,
#     zero = false,
# )

# figpath = "C:/Users/neutr/OneDrive - Imperial College London/PhD/Code/PolaronPathIntegrals/examples/plots/"

# for i in ["pdf", "tex", "svg"]
    # savefig(replicated_figure_one, figpath * "FHIP_replicated_fig_1.$i")
    # savefig(replicated_figure_two, figpath * "FHIP_replicated_fig_2.$i")
    # savefig(replicated_figure_three, figpath * "FHIP_replicated_fig_3.$i")
#     savefig(extended_figure_one, figpath * "FHIP_extend_fig_1_short.$i")
#     savefig(extended_figure_two, figpath * "FHIP_extend_fig_2_short.$i")
#     savefig(extended_figure_three, figpath * "FHIP_extend_fig_3_short.$i")
# end
