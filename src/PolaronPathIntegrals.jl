# PolaronPathIntegrals.jl

module PolaronPathIntegrals

using Optim
using QuadGK
using SpecialFunctions
using BigCombinatorics

include("frohlich.jl")
include("feynmantheory.jl")
include("osaka.jl")
include("hellwarththeory.jl")
include("mobility.jl")
include("optical_absorption.jl")
include("polaron.jl")

export frohlich_α
export feynman_variation
export feynman_free_energy
export osaka_free_energy
export singlemode_variation
export polaron_mobility
export polaron_mobility_zero
export optical_absorption
export optical_absorption_zero
export make_polaron

include("../src/polaronmakie/PolaronMakie.jl")
for name in names(PolaronMakie)
    @eval import .PolaronMakie: $(name)
    @eval export $(name)
end

# Physical constants
const ħ = 1.05457162825e-34; # kg m^2 s^{-1}
const eV = 1.602176487e-19; # kg m^2 s^{-2}
const m_e = 9.10938188e-31; # kg
const k_B = 1.3806504e-23; # kg m^2 K^{-1} s^2
const ϵ_0 = 8.854e-12 # C^2 N^{-1} m^{-2}
const c = 2.99792458e8 # m s^{-1}

struct Polaron
    α      # Frohlich alpha
    v      # Variational parameter
    w      # Variational parameter
    β      # Thermodynamic beta
    F      # Free energy
    μ      # Mobility
    Γ      # Absorption coefficient
end

function Base.show(io::IO, x::Polaron)
    print(io, "---------------------- \n Polaron Information: \n----------------------\n", "α = ", round(x.α, digits = 3), "\nv = ", round(x.v, digits = 3), "\nw = ", round(x.w, digits = 3), "\nβ = ", round(x.β, digits = 3), "\nF = ", round(x.F, digits = 3), "\nμ = ", round(Float64(x.μ), digits = 3), "\nΓ = ", round(Float64(x.Γ), digits = 3))
end

end # module
