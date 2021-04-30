#X(u)_function.jl

"""
Implementation of FHIP's susceptibility formula χ(Ω) for the Optical Absorption of Polarons and uses parts of Devresse's, Sitter's & Goovaerts paper that expands upon calculations in the FHIP paper (particularly for ℜχ(Ω)).

See FHIP 1962, equations (35a), (44), (47a, b, c):
https://link.aps.org/doi/10.1103/PhysRev.127.1004

and Devreese's et al. 1971 paper, equations (15), (16), and Appendix A:
https://link.aps.org/doi/10.1103/PhysRevB.5.2367.
"""

using QuadGK # Have to include this as this script somehow doesn't pick it up in the main PolaronMakie.jl module.

"""
----------------------------------------------------------------------
Simple implementation for finite temperatures. Useful for most values of arguments.
----------------------------------------------------------------------
"""

"""
ℜχ(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the real part of χ(Ω) (equation (35a) in FHIP) for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. β is the reduced thermodynamical beta = ħω/kT with k being Boltzman's constant. This attempts to just solve the integral directly using adaptive gaussian quadrature methods.
"""
function ℜχ(Ω, β, α, v, w)
    integrand(x) = (1 - cos(Ω * x)) * ℑS(x, 0.0, v, w, β, α)
    result = QuadGK.quadgk(x -> integrand(x), 0.0, Inf; atol = 1e-3)
    println(Ω, " ", result)
    return result
end

"""
ℑχ(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the imaginary part of χ(Ω) (equations (35a) and (44) in FHIP) for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. β is the reduced thermodynamical beta = ħω/kT with k being Boltzman's constant. This attempts to just solve the integral directly using adaptive gaussian quadrature methods.
"""
function ℑχ(Ω, β, α, v, w)
    integrand(x) = sin(Ω * x) * ℑS(x, 0.0, v, w, β, α)
    result = QuadGK.quadgk(x -> integrand(x), 0.0, Inf; atol = 1e-3)
    println(Ω, " ", result)
    return result
end

Ω_range = 0.01:0.1:10.01
ReX = [ℜχ(Ω, 0.4, 5, 4, 2.8) for Ω in Ω_range]
ImX = [ℑχ(Ω, 0.4, α, v, w) for Ω in Ω_range]
using Plots
plotly()
Plots.PlotlyBackend()
p = plot(Ω_range, [i[1] for i in ReX], ribbon = [i[2] for i in ReX], label = "ℜχ", xlabel = "Ω", ylabel = "χ")
plot!(p, Ω_range, [i[1] for i in ImX], ribbon = [i[2] for i in ImX], label = "ℑχ")
display(p)


# Create a complex number from the real and imaginary parts of χ(Ω).

function χ(Ω, β, α, v, w)

    R = (v^2 - w^2) / (w^2 * v)
    a_squared = β^2 / 4 + R * β * coth(β * v / 2)
    b = R * β / sinh(β * v / 2)

    D(x, y) = w^2 * (a_squared - β^2 / 4 - b * cos(v * x) * cosh(v * (y - β / 2)) + x^2 + y * (β - y)) / (β * v^2) + 1im * (w^2 * (b * sin(v * x) * sinh(v * (y - β / 2)) + 2 * x * (y - β / 2)) / (β * v^2))

    θ(x, y) = angle(D(x, y))
    r_squared(x, y) = abs2(D(x, y))

    S(x, y) = 2 * α * (cos(3 * θ(x, y) / 2) * cos(x) * cosh(y - β / 2) - sin(3 * θ(x, y) / 2) * sin(x) * sinh(y - β / 2)) / (3 * sqrt(π) * sinh(β / 2) * r_squared(x, y)^(3 / 4)) -2im * α * (cos(3 * θ(x, y) / 2) * sin(x) * sinh(y - β / 2) + sin(3 * θ(x, y) / 2) * cos(x) * cosh(y - β / 2)) / (3 * sqrt(π) * sinh(β / 2) * r_squared(x, y)^(3 / 4))

    integrand(x) = (1 - exp(1im * Ω * x)) * imag(S(x, 0.0))
    result = QuadGK.quadgk(x -> integrand(x), 0.0, Inf; atol = 1e-3)[1]
    println(Ω, " ", result)
    return result
end

function Γ(Ω, β, α, v, w)
    z = 1im * (Ω - χ(Ω, β, α, v, w) / Ω)
    abs(1 / z)
end
