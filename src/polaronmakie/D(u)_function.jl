# D(u)_function.jl

"""
Implementation of FHIP's susceptibility formula D(u) for the Optical Absorption of Polarons.

See FHIP 1962, equation (35c):
https://link.aps.org/doi/10.1103/PhysRev.127.1004.
"""

"""
ℜD(x::Float64, y::Float64, v::Float64, w::Float64, β::Float64)

    Calculate the real part of D(u) (equation (35c) in FHIP) for a complex argument u = x + iy. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. β is the reduced thermodynamical beta = ħω/kT with k being Boltzman's constant.
"""
function ℜD(x, y, v, w, β)
    x = BigFloat(x)
    y = BigFloat(y)
    β = BigFloat(β)
    v = BigFloat(v)
    w = BigFloat(w)

    R = (v^2 - w^2) / (w^2 * v)
    a_squared = β^2 / 4 + R * β * coth(β * v / 2)
    b = R * β / sinh(β * v / 2)

    w^2 * (a_squared - β^2 / 4 - b * cos(v * x) * cosh(v * (y - β / 2)) + x^2 + y * (β - y)) / (β * v^2)
end

"""
ℑD(x::Float64, y::Float64, v::Float64, w::Float64, β::Float64)

    Calculate the imaginary part of D(u) (equation (35c) in FHIP) for a complex argument u = x + iy. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. β is the reduced thermodynamical beta = ħω/kT with k being Boltzman's constant.
"""
function ℑD(x, y, v, w, β)
    x = BigFloat(x)
    y = BigFloat(y)
    β = BigFloat(β)
    v = BigFloat(v)
    w = BigFloat(w)

    R = (v^2 - w^2) / (w^2 * v)
    a_squared = β^2 / 4 + R * β * coth(β * v / 2)
    b = R * β / sinh(β * v / 2)

    w^2 * (b * sin(v * x) * sinh(v * (y - β / 2)) + 2 * x * (y - β / 2)) / (β * v^2)
end

# Create a complex number from the real and imaginary parts of D(u).

D(x, y, v, w, β) = ℜD(x, y, v, w, β) + 1im * ℑD(x, y, v, w, β)
