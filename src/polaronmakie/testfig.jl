using GLMakie
using PolaronPathIntegrals

# Physical constants
const ħ = 1.05457162825e-34; # Reduced Planck's constant (kg m^2 s^{-1})
const eV = 1.602176487e-19; # Electron Volt (kg m^2 s^{-2})
const m_e = 9.10938188e-31; # Electron Mass (kg)
const k_B = 1.3806504e-23; # Boltzmann's constant (kg m^2 K^{-1} s^2)
const ϵ_0 = 8.854e-12 # Dielectric constant (C^2 N^{-1} m^{-2})
const c = 2.99792458e8 # Speed of light (m s^{-1})

fig = Figure(resolution = (1800, 900), backgroundcolor = RGBf0(0.98, 0.98, 0.98))

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

set_close_to!(lsgrid.sliders[1], 24.1)
set_close_to!(lsgrid.sliders[2], 4.5)
set_close_to!(lsgrid.sliders[3], 0.12)
set_close_to!(lsgrid.sliders[4], 2.25)

fig[3, 1:5] = lsgrid.layout


Ω_slider = Slider(fig[1:2, 1], range = 0.001:0.1:10.001, startvalue = 0.001, horizontal = false, tellwidth = true, height = nothing, width = Auto())
freq_sliders = [s.value for s in vcat(lsgrid.sliders , [Ω_slider])]
Ω_slider_label = Label(fig[1:2, 0], lift(s -> "Ω: $(round(s[], digits = 2)) Hz", Ω_slider.value), rotation = pi/2)

polaron_temp = lift(freq_sliders...) do slvalues...
    make_polaron(slvalues[2], slvalues[1], slvalues[4] * 1e12, slvalues[3]; temp = 100:10:400, efield_freq = slvalues[5], verbose = false)
end

T = lift(p -> p.T, polaron_temp)
v = lift(p -> Vector{Float64}(p.v), polaron_temp)
w = lift(p -> Vector{Float64}(p.w), polaron_temp)
κ = lift(p -> Vector{Float64}(p.κ), polaron_temp)
M = lift(p -> Vector{Float64}(p.M), polaron_temp)
F = lift(p -> Vector{Float64}(p.F), polaron_temp)
μ_T = lift(p -> vcat.(p.μ...)[1], polaron_temp)
Γ_T = lift(p -> vcat.(p.Γ...)[1], polaron_temp)

T_slider = Slider(fig[1:2, 7], range = 100:1:400, startvalue = 300, horizontal = false, tellwidth = true, height = nothing, width = Auto())
T_slider_label = Label(fig[1:2, 8], lift(s -> "T: $(s[]) K", T_slider.value), rotation = -pi/2)

temp_sliders = [s.value for s in vcat(lsgrid.sliders , [T_slider])]

polaron_freq = lift(temp_sliders...) do slvalues...
    make_polaron(slvalues[2], slvalues[1], slvalues[4] * 1e12, slvalues[3]; temp = slvalues[5], efield_freq = 0.001:0.5:10.001, verbose = false)
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

# μ_Ω = Vector{Any}()
# Γ_Ω = Vector{Any}()
# for i in 1:Int(floor(length(T[])/(N - 1))):length(T[])
#     push!(μ_Ω, lift(p -> Vector{Float64}(p.μ[i]), polaron))
#     push!(Γ_Ω, lift(p -> Vector{Float64}(p.Γ[i]), polaron))
# end
# μ_T = Vector{Any}()
# Γ_T = Vector{Any}()
# for i in 1:Int(floor(length(Ω[])/(N - 1))):length(Ω[])
#     push!(μ_T, lift(p -> vcat.(p.μ...)[i], polaron))
#     push!(Γ_T, lift(p -> vcat.(p.Γ...)[i], polaron))
# end
# for (μ, Γ) in zip(μ_T, Γ_T)
#     scatterlines!(mobility_axis, T, μ, colormap = :viridis)
#     scatterlines!(absorption_axis, T, Γ, color = :blue)
# end