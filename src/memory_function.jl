# memory_function.jl

"""
χ(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) at finite temperature for a given frequency Ω. β is the thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""

function χ(Ω, β, α, v, w)
    R = (v^2 - w^2) / (w^2 * v)
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β * v / 2) + x^2 / β - 1im * (R * sin(v * x) + x))
    S(x) = 2 * α / (3 * √π) * (exp(1im * x) + 2 * cos(x) / (exp(β) - 1)) / (D(x))^(3 / 2)
    integrand(x) = (1 - exp(-1im * Ω * x)) * imag(S(x))
    result = QuadGK.quadgk(x -> integrand(x), 0.0, Inf)[1]
    return result
end

"""
χ(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) at zero temperature for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""

function χ(Ω, α, v, w)
    if Ω < 1
        return 0
    else
        R = (v^2 - w^2) / (w^2 * v)
        D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) + (x / 20)^2 - 1im * (R * sin(v * x) + x))
        S(x) = 2 * α / (3 * √π) * exp(1im * x) / (D(x))^(3 / 2)
        integrand(x) = (1 - exp(-1im * Ω * x)) * imag(S(x))
        QuadGK.quadgk(x -> integrand(x), 0.0, Inf)[1]
    end
end

"""
function χ(Ω, β, α, v, w)

    R = (v^2 - w^2) / (w^2 * v)
    a_squared = β^2 / 4 + R * β * coth(β * v / 2)
    b = R * β / sinh(β * v / 2)

    D(x, y) = w^2 * (a_squared - β^2 / 4 - b * cos(v * x) * cosh(v * (y - β / 2)) + x^2 + y * (β - y)) / (β * v^2) + 1im * (w^2 * (b * sin(v * x) * sinh(v * (y - β / 2)) + 2 * x * (y - β / 2)) / (β * v^2))

    θ(x, y) = angle(D(x, y))
    r_squared(x, y) = abs2(D(x, y))

    S(x, y) = 2 * α * (cos(3 * θ(x, y) / 2) * cos(x) * cosh(y - β / 2) - sin(3 * θ(x, y) / 2) * sin(x) * sinh(y - β / 2)) / (3 * sqrt(π) * sinh(β / 2) * r_squared(x, y)^(3 / 4)) -2im * α * (cos(3 * θ(x, y) / 2) * sin(x) * sinh(y - β / 2) + sin(3 * θ(x, y) / 2) * cos(x) * cosh(y - β / 2)) / (3 * sqrt(π) * sinh(β / 2) * r_squared(x, y)^(3 / 4))

    integrand(x) = (1 - exp(-1im * Ω * x)) * imag(S(x, 0.0))
    result = QuadGK.quadgk(x -> integrand(x), 0.0, Inf; atol = 1e-3)[1]
    # println(Ω, " ", result)
    return result
end
"""
