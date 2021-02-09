# mobility.jl

"""
----------------------------------------------------------------------
Finite temperature implementation for ℑχ using BesselK functions.
----------------------------------------------------------------------

Details come from Appendix A of Devreese's et al. paper, although they take the limit of β -> ∞, whereas here we do not. Likewise, we provide proper treatment of any expansions without immediate approximatons (i.e. in Appendix A they ignore the even part of the cos^n(vx) expansion without any apparent justification).
"""

"""
ℑχ(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the imaginary part of χ(Ω) in a zero temperature approximation (equation (16) in Devreese's et al.) for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function ℑχ(Ω, α, v, w)

    if Ω < 1
        return 0.0
    else
        R = (v^2 - w^2) / (w^2 * v)
        coefficient = 2 / 3 * α * (v / w)^3

        total_sum = 0.0
        next_current = 0.0
        n = 0

        while next_current >= total_sum * 1e-8
            current =
                -coefficient *
                sqrt(π) *
                gamma(n + 1) *
                (-2 * R)^n *
                abs(Ω - 1 - n * v)^(n + 1 / 2) *
                exp(-R * abs(Ω - 1 - n * v)) *
                (1 + sign(Ω - 1 - n * v)) / (
                    gamma(-n - 1 / 2) *
                    BigCombinatorics.doublefactorial(2 * n + 1)
                )
            next_current = current
            total_sum += next_current
            n += 1
        end
        return total_sum
    end
end

"""
ℑχ(Ω::Float64, α::Float64, v::Float64, w::Float64, β::Float64)

    Calculate the imaginary part of χ(Ω) at finite temperature using an infinite expansion of BesselK functions (see Appendix A in Devreese's et al.), for a given frequency Ω. β is the thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""

function ℑχ(Ω, β, α, v, w)

    # Set arguments to BigFloat precision. Without this the calculations break down due to large values of hyperbolic functions.
    Ω = BigFloat(Ω)
    β = BigFloat(β)
    α = BigFloat(α)
    v = BigFloat(v)
    w = BigFloat(w)

    # Initialise constants.
    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    # Coefficient independent of summation variables n & k.
    coefficient =
        2 * α / (3 * sqrt(π)) * (v / w)^3 * β^(3 / 2) * sinh(Ω * β / 2) /
        sinh(β / 2)

    # Initialise total value of double summation as a BigFloat.
    total_sum = BigFloat(0.0)
    next_current = BigFloat(0.0)
    n = 0

    while next_current >= total_sum * 1e-3 # Infinite summation from the Binomial expansion of (x^2 + a^2 - b * cos(v * x))^(-3/2). |b * cos(v * x) / (x^2 + a^2)| < 1 for v > 0 and β > 0. Here just limit summation to when next term adds negiglible amount to total sum.

        current = BigFloat(0.0)

        # Coefficient that depends on n summation variable.
        n_coefficient = -2 * SpecialFunctions.beta(-1 / 2 - n, n + 1)^(-1) * (-1)^n * b^n * (4 * a)^(-n - 1) * sqrt(π) / SpecialFunctions.gamma(n + 3 / 2)

        # Summation over expansion of (cos(v * x))^n splits into sums over odd n and even n.

        # Finite sum over k for even n from (cos(v * x))^n expansion.
        for k = -1:floor((n - 1) / 2)

            # Coefficient that depends on k summation variable.
            k_coefficient_inverse = vcat(
                repeat([(n + 1) * SpecialFunctions.beta(n - k + 1, k + 1)], 4),
                repeat([(n + 1) * SpecialFunctions.beta(n / 2 + 1, n / 2 + 1)], 2,),
            )

            # Arguments of resultant cosines obtained from cos^n(vx) * cos(x) * cos(Ωx).
            y = [
                (Ω + 1 + v * (n - 2 * k)),
                (Ω + 1 - v * (n - 2 * k)),
                (Ω - 1 + v * (n - 2 * k)),
                (Ω - 1 - v * (n - 2 * k)),
                (Ω + 1),
                (Ω - 1),
            ]

            for (x, k_coeff) in zip(y, k_coefficient_inverse)  # Iterate over cosine arguments for brevity.

                """
                Integral over cos(yx) / (x^2 + a^2)^(n + 3/2) can be identified with BesselK functions.
                """

                # Sometimes SpecialFunctions.besselk overflow, so catch error if it does and instead use ArbNumerics.besselk.
                err = Nothing
                bessel = try
                    SpecialFunctions.besselk(Int(n + 1), abs(Float64(x)) * Float64(a))
                catch err
                    err
                end

                # If SpecialFunctions.besselk overflows or loses accuracy (i.e. equates to 0.0 which shouldn't happen in this range) use arbitrary precision.
                if err == SpecialFunctions.AmosException(2) || bessel == 0.0

                    setextrabits(0) # Do not use extra bits for evals.
                    p = Int64(precision(coefficient)) # Start at BigFloat prec.
                    bessel = ArbNumerics.besselk(n + 1, abs(ArbFloat(x * a, bits = p)))

                    # Increment precision until overflows / loss of accuracy ceases.
                    while bessel == NaN || bessel == 0.0
                        bessel = ArbNumerics.besselk(n + 1, abs(ArbFloat(x * a, bits = p)))
                        p += 1
                    end
                end

                # Calculate nth term in expansion in terms of BesselK functions.
                current += bessel * abs(x)^(n + 1) / k_coeff
            end
        end

        # Remember nth term for convergence comparison and add to total sum.
        next_current = n_coefficient * current
        total_sum += n_coefficient * current

        # Go to next n.
        n += 1
    end

    # Return final value obtained from double summation.
    return coefficient * total_sum
end

"""
----------------------------------------------------------------------
The Moblity of the Polaron.
----------------------------------------------------------------------
"""

"""
polaron_mobility(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the moblity μ(Ω) of the polaron at zero-temperatures (equation (1) in Hellwarth 1999) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function polaron_mobility(Ω, α, v, w)
    Ω / ℑχ(Ω, α, v, w)
end

"""
polaron_mobility(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the moblity μ(Ω) of the polaron at finite temperatues (equation (1) in Hellwarth 1999) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function polaron_mobility(Ω, α, v, w, β)
    Ω / ℑχ(Ω, α, v, w, β)
end
