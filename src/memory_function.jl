# memory_function.jl

"""
χ(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) of the polaron at finite temperatures (equation (35a) in FHIP 1962) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function χ(Ω, β, α, v, w)

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)

    # FHIP1962, page 1009, eqn (35c).
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β * v / 2) + x^2 / β - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36).
    S(x) = 2 * α / (3 * √π) * (exp(1im * x) + 2 * cos(x) / (exp(β) - 1)) / (D(x))^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(x) = (1 - exp(-1im * Ω * x)) * imag(S(x))
    QuadGK.quadgk(x -> integrand(x), 0.0, Inf)[1]
end

"""
χ(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function χ(Ω) of the polaron at zero-temperatures (equation (35a) in FHIP 1962) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function χ(Ω, α, v, w)

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)

    # FHIP1962, page 1009, eqn (35c) with β → ∞.
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36) with β → ∞.
    S(x) = 2 * α / (3 * √π) * exp(1im * x) / (D(x))^(3 / 2)

    # FHIP1962, page 1009, eqn (35a). Set upper limit < ∞ but very large so cos argument finite.
    integrand(x) = (1 - exp(-1im * Ω * x)) * imag(S(x))
    QuadGK.quadgk(x -> integrand(x), 0.0, 1e200)[1]
end

"""
χ_dc(β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the memory function lim(Ω → 0){χ(Ω) / Ω} of the polaron at finite temperatures (equation (35a) in FHIP 1962) at zero frequency. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function χ_dc(β, α, v, w)

    # FHIP1962, page 1011, eqn (47c).
    R = (v^2 - w^2) / (w^2 * v)

    # FHIP1962, page 1009, eqn (35c).
    D(x) = w^2 / v^2 * (R * (1 - cos(v * x)) * coth(β * v / 2) + x^2 / β - 1im * (R * sin(v * x) + x))

    # FHIP1962, page 1009, eqn (36).
    S(x) = 2 * α / (3 * √π) * (exp(1im * x) + 2 * cos(x) / (exp(β) - 1)) / (D(x))^(3 / 2)

    # Set frequency small enough to mimic Ω = 0 without generating numerical instabilities in integral.
    Ω = 1e-200

    # FHIP1962, page 1009, eqn (35a). Readily divided by Ω to get sinc function in integrand rather than sine which would give 0 all the time.
    integrand(x) = (1 - exp(-1im * Ω * x)) * imag(S(x)) / Ω
    QuadGK.quadgk(x -> integrand(x), 0.0, Inf)[1]
end

# for i in 1:length(α)
#     append!(s, [χ(i, βr[j], α[j], 3.36600719617118, 2.6419291934269666) for i in 0.1:0.1:3.5])
# end