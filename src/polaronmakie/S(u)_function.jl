# S(u)_function.jl

"""
Implementation of FHIP's susceptibility formula S(u) for the Optical Absorption of Polarons.

See FHIP 1962, equation (36):
https://link.aps.org/doi/10.1103/PhysRev.127.1004.
"""

# Calculate the modulus-squared of the complex function D(u) from equation (35c) in FHIP.

mod_squared_D(x, y, v, w, β) = abs2(D(x, y, v, w, β))

# Calculate the complex argument of the complex function D(u) from equation (35c) in FHIP.

arg_D(x, y, v, w, β) = angle(D(x, y, v, w, β))

"""
ℜS(x::Float64, y::Float64, v::Float64, w::Float64, β::Float64, α::Float64)

    Calculate the real part of S(u) (equation (36) in FHIP) for a complex argument u = x + iy. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. β is the reduced thermodynamical beta = ħω/kT with k being Boltzman's constant.
"""
function ℜS(x, y, v, w, β, α)
    x = BigFloat(x)
    y = BigFloat(y)
    β = BigFloat(β)
    α = BigFloat(α)
    v = BigFloat(v)
    w = BigFloat(w)

    θ = angle(D(x, y, v, w, β))
    r_squared = abs2(D(x, y, v, w, β))

    2 * α * (cos(3 * θ / 2) * cos(x) * cosh(y - β / 2) - sin(3 * θ / 2) * sin(x) * sinh(y - β / 2)) / (3 * sqrt(π) * sinh(β / 2) * r_squared^(3 / 4))
end

"""
ℑS(x::Float64, y::Float64, v::Float64, w::Float64, β::Float64, α::Float64)

    Calculate the imaginary part of S(u) (equation (36) in FHIP) for a complex argument u = x + iy. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. β is the reduced thermodynamical beta = ħω/kT with k being Boltzman's constant.
"""
function ℑS(x, y, v, w, β, α)
    x = BigFloat(x)
    y = BigFloat(y)
    β = BigFloat(β)
    α = BigFloat(α)
    v = BigFloat(v)
    w = BigFloat(w)

    θ = angle(D(x, y, v, w, β))
    r_squared = abs2(D(x, y, v, w, β))

    -2 * α * (cos(3 * θ / 2) * sin(x) * sinh(y - β / 2) + sin(3 * θ / 2) * cos(x) * cosh(y - β / 2)) / (3 * sqrt(π) * sinh(β / 2) * r_squared^(3 / 4))
end

# Create a complex number from the real and imaginary parts of S(u).

S(x, y, v, w, β, α) = ℜS(x, y, v, w, β, α) + 1im * ℑS(x, y, v, w, β, α)
