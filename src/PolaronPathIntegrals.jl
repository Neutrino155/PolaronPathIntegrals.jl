# PolaronPathIntegrals.jl

module PolaronPathIntegrals

using Optim
using QuadGK

include("frohlich.jl")
include("feynmantheory.jl")
include("osaka.jl")
include("hellwarththeory.jl")
include("mobility.jl")
include("optical_absorption.jl")


export frohlich_α
export feynman_variation
export feynman_free_energy
export osaka_free_energy
export singlemode_variation
export polaron_mobility
export polaron_mobility_zero
export optical_absorption
export optical_absorption_zero

include("../src/polaronmakie/PolaronMakie.jl")
for name in names(PolaronMakie)
    @eval import .PolaronMakie: $(name)
    @eval export $(name)
end

# Physical constants
const ħ = 1.05457162825e-34; # kg m^2 s^{-1
const eV = 1.05457162825e-34; # kg m^2 s^{-2}
const m_e = 9.10938188e-31; # kg
const k_B = 1.3806504e-23; # kg m^2 K^{-1} s^2
const ϵ_0 = 8.854e-12 # C^2 N^{-1} m^{-2}
const c = 2.99792458e8 # m s^{-1}

struct Polaron
    α::Float64      # Frohlich alpha
    v::Float64      # Variational parameter
    w::Float64      # Variational parameter
    F::Float64      # Free energy
    μ::Float64      # Mobility 
end

function Polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff; temp = 0.0, field_freq = 0.0)
    ω = 2 * π * phonon_freq
    α = frohlich_α(ϵ_optic, ϵ_static, phonon_freq, m_eff)
    if temp = 0.0
        v, w = feynman_variation(α; v = 7.2, w = 6.5)
        F = feynman_free_energy(v, w, α)
        μ = ω * m_eff * polaron_mobility_zero(e_freq, α, v, w, N = 10) / eV
    else
        β = ħ * ω / (k_B * temp)
        v, w = single_variation(α, β; v = 7.2, w = 6.5)
        F = osaka_free_energy(v, w, β, α)
        μ = ω * m_eff * polaron_mobility(e_freq, β, α, v, w, N = 10) / eV
    end
    new(α, v, w, β, F, μ)
end

function Base.show(io::IO, quantities::Polaron)
    print(io, "$α, $v, $w, $β, $F, $μ")
end
