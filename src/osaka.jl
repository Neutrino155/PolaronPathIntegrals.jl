#osaka.jl

"""
    Implementation of Osaka's extended treatment of Feynman's variational technique applied to the Polaron model; generalising it from the case at O^{∘}K to the case at finite temperature.

    See Osaka 1959:
    Progress of Theoretical Physics, Vol. 22, No.3, September 1959

    And Hellwarth 1999:
    https://doi.org/10.1103/PhysRevB.60.299

"""

# Equation 62b in Hellwarth. Equation 20 in Osaka.
A(v, w, β) = 3 / β * (log(v / w) - 1 / 2 * log(2 * π * β) - log(sinh(v * β / 2) / sinh(w * β / 2)))

# Equation 62d in Hellwarth.
Y(x, v, β) = 1 / (1 - exp(-v * β)) * (1 + exp(-v * β) - exp(-v * x) - exp(v * (x - β)))

# Osaka's version: Within Equation 16.
# Y(x, v, β) = 1 - exp(-v * x) + (1 - coth(v * β / 2)) * (cosh(v * x) - 1)

# Integrand of Equation 62c in Hellwarth.
B_integrand(x, v, w, β) =  (exp(β - x) + exp(x)) / sqrt(w^2 * x * (1 - x / β) + Y(x, v, β) * (v^2 - w^2) / v)

# Osaka's version: Within Equation 16.
# B_integrand(x, v, w, β) =  (exp(-x)) / sqrt(abs(w^2 * x * (1 - x / β) + Y(x, v, β) * (v^2 - w^2) / v))

# Equation 62c in Hellwarth.
B(v, w, β, α) = α * v / (sqrt(π) * (exp(β) - 1)) * QuadGK.quadgk(x -> B_integrand(x, v, w, β), 0, β / 2)[1]

# Osaka's version: Equation 16.
# B(v, w, β, α) = α * v * exp(β) / (sqrt(π) * (exp(β) - 1)) * QuadGK.quadgk(x -> B_integrand(x, v, w, β), 0, β)[1]

# Equation 62e in Hellwarth. Equation 17 in Osaka.
C(v, w, β) = 3 / 4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β))

# Equation 62a in Hellwarth. In paragraph below Equation 22 in Osaka; has extra 1/β due to different definition of A, B & C.
osaka_free_energy(v, w, β, α) = -(A(v, w, β) + B(v, w, β, α) + C(v, w, β))
