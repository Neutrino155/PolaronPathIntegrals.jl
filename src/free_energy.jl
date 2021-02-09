# free_energy.jl

"""
Implementation of Feynman's original variational technique applied to the Polaron model.

See Feynman 1955:
http://dx.doi.org/10.1103/PhysRev.97.660
"""

# Equation 31: The <|X(t) - X(s)|^{-1}> * exp(-|t-w|) effective action.
A_integrand(x, v, w) = (w^2 * x + (v^2 - w^2) / v * (1 - exp(-v * x)))^(-0.5) * exp(-x)

A(v, w, α) = π^(-0.5) * α * v * QuadGK.quadgk(x -> A_integrand(x, v, w), 0, Inf)[1]

# Equation 33: Lowest Free energy E = -B - A where B = -3/(4v)*(v-w)^2.
free_energy(v, w, α) = (3 / (4 * v)) * (v - w)^2 - A(v, w, α)

"""
    Implementation of Osaka's extended treatment of Feynman's variational technique applied to the Polaron model; generalising it from the case at O^{∘}K to the case at finite temperature.

    See Osaka 1959:
    Progress of Theoretical Physics, Vol. 22, No.3, September 1959

    And Hellwarth 1999:
    https://doi.org/10.1103/PhysRevB.60.299

"""

# Equation 62d in Hellwarth.
Y(x, v, β) = 1 / (1 - exp(-v * β)) * (1 + exp(-v * β) - exp(-v * x) - exp(v * (x - β)))

# Integrand of Equation 62c in Hellwarth.
A_integrand(x, v, w, β) =  (exp(β - x) + exp(x)) / sqrt(w^2 * x * (1 - x / β) + Y(x, v, β) * (v^2 - w^2) / v)

# Equation 62c in Hellwarth.
err = false
A(v, w, α, β) = try
    α * v / (sqrt(π) * (exp(β) - 1)) * QuadGK.quadgk(x -> A_integrand(x, v, w, β), 0, β / 2)[1]
catch
    err = true
    @show(A(v, w, α, β))
end

if err
    A(v, w, α, β) = α * v / (sqrt(π) * (exp(BigFloat(β)) - 1)) * QuadGK.quadgk(x -> A_integrand(x, v, w, BigFloat(β)), 0, BigFloat(β / 2))[1]
    @show(A(v, w, α, β))
end


# Equation 62b in Hellwarth. Equation 20 in Osaka.
B(v, w, β) = 3 / β * (log(v / w) - 1 / 2 * log(2 * π * β) - log(sinh(v * β / 2) / sinh(w * β / 2)))

# Equation 62e in Hellwarth. Equation 17 in Osaka.
C(v, w, β) = 3 / 4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β))

# Equation 62a in Hellwarth. In paragraph below Equation 22 in Osaka; has extra 1/β due to different definition of A, B & C.
free_energy(v, w, α, β) = -(A(v, w, α, β) + B(v, w, β) + C(v, w, β))
