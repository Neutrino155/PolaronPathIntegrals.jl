# optical_absorption.jl

"""
----------------------------------------------------------------------
Polaron absorption coefficient Γ(Ω).
----------------------------------------------------------------------
"""

"""
optical_absorption_ac(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the absorption coefficient Γ(Ω) for the polaron at finite temperatures (equation (11a) in Devreese's et al.) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function optical_absorption_ac(Ω, β, α, v, w)
    X = χ(Ω, β, α, v, w)
    Reχ = real(X)
    Imχ = imag(X)
    result = Ω * Imχ / (Ω^4 - 2 * Ω^2 * Reχ + Reχ^2 + Imχ^2)
    return result
end

"""
optical_absorption_ac(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the absorption coefficient Γ(Ω) for the polaron at zero temperatures (equation (11a) in Devreese's et al.) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function optical_absorption_ac(Ω, α, v, w)
    X = χ(Ω, α, v, w)
    Reχ = real(X)
    Imχ = imag(X)
    result = Ω * Imχ / (Ω^4 - 2 * Ω^2 * Reχ + Reχ^2 + Imχ^2)
    return result
end

"""
optical_absorption_dc(β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the absorption coefficient Γ(Ω) for the polaron at finite temperatures (equation (11a) in Devreese's et al.) at zero frequency. β is the thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function optical_absorption_dc(β, α, v, w)
    X = χ_dc(β, α, v, w)
    Reχ = real(X)
    Imχ = imag(X)
    result = Ω * Imχ / (Ω^4 - 2 * Ω^2 * Reχ + Reχ^2 + Imχ^2)
    return result
end