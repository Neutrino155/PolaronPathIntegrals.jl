#frohlich.jl

"""
Implementation of Frohlich's Hamiltonian and coupling parameter α applied to the Polaron model.

See Frohlich 1952:
https://www.jstor.org/stable/99166

"""

"""
frohlich_α(ϵ_optic::Float64, ϵ_static::Float64, freq::Float64, m_eff::Float64)

    Calculates the Frohlich alpha parameter, for a given dielectric constant, frequency (f) of phonon in Hertz, and effective mass (in units of the bare electron mass).

    See Equation 2.37.

"""
function frohlich_α(ϵ_optic::Float64, ϵ_static::Float64, freq::Float64, m_eff::Float64)

    ω = 2 * π * freq

    α = 0.5 / (4 * π * ϵ_0) *
       (1 / ϵ_optic - 1 / ϵ_static) *
       (eV^2 / (ħ * ω)) *
       sqrt(2 * m_e * m_eff * ω / ħ)

    return α
end
