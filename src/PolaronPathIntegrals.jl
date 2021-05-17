# PolaronPathIntegrals.jl

module PolaronPathIntegrals

using Optim
using QuadGK
using Plots
using AbstractPlotting
using GLMakie
using PyPlot

include("coupling.jl")
include("free_energy.jl")
include("variation.jl")
include("memory_function.jl")
include("mobility.jl")
include("optical_absorption.jl")
include("make_polaron.jl")
include("plot_polaron.jl")

export frohlich_α
export variation
export free_energy
export polaron_mobility
export optical_absorption
export make_polaron
export plot_polaron
export save_polaron_plots
export plot_polaron_interactive

# Individual Functions
# export χ
include("../src/polaronmakie/PolaronMakie.jl")
for name in names(PolaronMakie)
    @eval import .PolaronMakie: $(name)
    @eval export $(name)
end

# Physical constants
const ħ = 1.05457162825e-34; # Reduced Planck's constant (kg m^2 s^{-1})
const eV = 1.602176487e-19; # Electron Volt (kg m^2 s^{-2})
const m_e = 9.10938188e-31; # Electron Mass (kg)
const k_B = 1.3806504e-23; # Boltzmann's constant (kg m^2 K^{-1} s^2)
const ϵ_0 = 8.854e-12 # Dielectric constant (C^2 N^{-1} m^{-2})
const c = 2.99792458e8 # Speed of light (m s^{-1})
const amu = 1.66053906660e-27 # Atomic Mass Unit (kg)

struct Polaron
    α      # Frohlich alpha (unitless)
    T      # Temperature (K)
    β      # Reduced Thermodynamic beta (unitless)
    v      # Variational parameter (s^-1)
    w      # Variational parameter (s^-1)
    κ      # Fictitious spring constant (multiples of m_e) (kg / s^2)
    M      # Fictitious particle (multiples of m_e) (kg)
    F      # Free energy (meV)
    Ω      # Electric field frequencies (multiples of phonon frequency ω) (s^-1)
    μ      # Mobility (cm^2 / Vs)
    Γ      # Absorption coefficient (cm^-1)
end

# Broadcast Polaron data.
function Base.show(io::IO, x::Polaron)
    flush(stdout)
    print(io, "---------------------- \n Polaron Information: \n----------------------\n", "α = ", round(x.α, digits = 3), "\nT = ", round.(x.T, digits = 3), " K \nβ = ", round.(x.β, digits = 3), "\nv = ", round.(x.v, digits = 3), " s^-1\nw = ", round.(x.w, digits = 3), " s^-1\nκ = ", round.(x.κ, digits = 3), " kg/s^2\nM = ", round.(x.M, digits = 3), " kg\nF = ", round.(x.F, digits = 3), " meV\nΩ = ", round.(Float64.(x.Ω), digits = 3),  " s^-1\nμ = ", x.μ .|> y -> round.(Float64.(y), digits = 3), " cm^2/Vs\nΓ = ", x.Γ .|> y -> round.(Float64.(y), digits = 3), " cm^-1")
end

export Polaron

end # module
