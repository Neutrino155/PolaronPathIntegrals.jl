"""
----------------------------------------------------------------------
Finite temperature implementation for ℑχ using BesselK functions.
----------------------------------------------------------------------

Details come from Appendix A of Devreese's et al. paper, although they take the limit of β -> ∞, whereas here we do not. Likewise, we provide proper treatment of any expansions without immediate approximatons (i.e. in Appendix A they ignore the even part of the cos^n(vx) expansion without any apparent justification).
"""

"""
ℑχ(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64, N::Int)

    Calculate the imaginary part of χ(Ω) at finite temperature using an infinite expansion of BesselK functions (see Appendix A in Devreese's et al.), for a given frequency Ω. β is the thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. N is the upper limit of a sum that is analytically infinite, however most of the information is encapuslated by N < 50.

    Function can either implement the BesselK functions to FLoat64 precision using the SpecialFunctions.jl package (which eventually Overflows), or to arbitrary precision with the ArbNumerics.jl package. Change comments where appropriate. Set the precision of the ArbReal type with 'setworkingprecision(ArbReal, bits=128)'. This needs to be the same precision as BigFloat, set with 'setprecision(BigFloat, 128)'.
"""
function ℑχ(Ω, β, α, v, w, N = 10) # larger N increases accuracy.

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
        2 * α / (3 * sqrt(π)) * (v / w)^3 * β^(3 / 2) * sinh(β * Ω / 2) / sinh(β / 2)

    # Initialise total value of double summation as a BigFloat.

    total_sum = BigFloat(0.0)

    for n = 0:Int(N)   # Infinite summation from the Binomial expansion of (x^2 + a^2 - b * cos(v * x))^(-3/2). |b * cos(v * x) / (x^2 + a^2)| < 1 for v > 0 and β > 0.

        # Coefficient that depends on n summation variable.

        n_coefficient =
            -2 * SpecialFunctions.beta(-1 / 2 - n, n + 1)^(-1) * (-1)^n * b^n * (4 * a)^(-n - 1) * sqrt(π) /
            SpecialFunctions.gamma(n + 3 / 2)

        # Summation over expansion of (cos(v * x))^n splits into sums over odd n and even n.

        if isodd(n) # For odd n.

            # Finite sum over k for odd n from (cos(v * x))^n expansion.

            for k = -1:Int((n - 1) / 2)

                # Coefficient that depends on k summation variable.

                k_coefficient_inverse = ((n + 1) * SpecialFunctions.beta(n - k + 1, k + 1))

                # Arguments of resultant cosines obtained from cos^n(vx) * cos(x) * cos(Ωx).

                y = [
                    (Ω + 1 + v * (n - 2 * k)),
                    (Ω + 1 - v * (n - 2 * k)),
                    (Ω - 1 + v * (n - 2 * k)),
                    (Ω - 1 - v * (n - 2 * k)),
                ]

                for x in y # Iterate over cosine arguments for brevity.

                    """
                    For low-ish values of β (~ < 40) and Ω (~ < 50). Integral over cos(yx) / (x^2 + a^2)^(n + 3/2) can be identified with BesselK functions.
                    """

                    total_sum +=
                        coefficient *
                        n_coefficient *
                        besselk(Int(n + 1), abs(Float64(x)) * Float64(a)) *
                        abs(x)^(n + 1) / k_coefficient_inverse

                    """
                    For larger values of β (~ > 40) and Ω (~ > 50). This is an arbitrary precision implementation of BesselK functions by ArbNumerics.jl that is not available in the SpecialFunctions.jl package. Comment above and uncomment below for higher precision calculations.
                    """

                    # Note: May be type issues somewhere resulting in Any types and slowing down code.

                    # total_sum += coefficient * n_coefficient * k_coefficient * besselk(ArbReal(n + 1), ArbReal(abs(x) * a)) * abs(x)^(n + 1)
                end
            end
        end

        if iseven(n) # For even n.

            # Finite sum over k for even n from (cos(v * x))^n expansion.

            for k = -1:Int(n / 2 - 1)

                # Coefficient that depends on k summation variable.

                k_coefficient_inverse = ((n + 1) * SpecialFunctions.beta(n - k + 1, k + 1))

                # Arguments of resultant cosines obtained from cos^n(vx) * cos(x) * cos(Ωx).

                y = [
                    (Ω + 1 + v * (n - 2 * k)),
                    (Ω + 1 - v * (n - 2 * k)),
                    (Ω - 1 + v * (n - 2 * k)),
                    (Ω - 1 - v * (n - 2 * k)),
                ]

                for x in y # Iterate over cosine arguments for brevity.

                    """
                    For low-ish values of β (~ < 40) and Ω (~ < 50). Integral over cos(yx) / (x^2 + a^2)^(n + 3/2) can be identified with BesselK functions.
                    """

                    total_sum +=
                        coefficient *
                        n_coefficient *
                        besselk(Int(n + 1), abs(Float64(x)) * Float64(a)) *
                        abs(x)^(n + 1) / k_coefficient_inverse

                    """
                    For larger values of β (~ > 40) and Ω (~ > 50). This is an arbitrary precision implementation of BesselK functions by ArbNumerics.jl that is not available in the SpecialFunctions.jl package. Comment above and uncomment below for higher precision calculations.
                    """

                    # Note: May be type issues somewhere resulting in Any types and slowing down code.

                    # total_sum += coefficient * n_coefficient * k_coefficient * besselk(ArbReal(n + 1), ArbReal(abs(x) * a)) * abs(x)^(n + 1)
                end

                # Coefficient that depends on k summation variable. This term is an extra term that arises for the even part of the cos^n(vx) expansion.

                k_even_coefficient_inverse = ((n + 1) * SpecialFunctions.beta(n / 2 + 1, n / 2 + 1))

                # Arguments of resultant cosines obtained from just cos(x) * cos(Ωx).

                y_even = [(Ω + 1), (Ω - 1)]

                for x in y_even # Iterate over cosine arguments for brevity.

                    """
                    For low-ish values of β (~ < 40) and Ω (~ < 50). Integral over cos(y_even x) / (x^2 + a^2)^(n + 3/2) can be identified with BesselK functions.
                    """

                    total_sum +=
                        coefficient *
                        n_coefficient *
                        besselk(Int(n + 1), abs(Float64(x)) * Float64(a)) *
                        abs(x)^(n + 1) / k_even_coefficient_inverse

                    """
                    For larger values of β (~ > 40) and Ω (~ > 50). This is an arbitrary precision implementation of BesselK functions by ArbNumerics.jl that is not available in the SpecialFunctions.jl package. Comment above and uncomment below for higher precision calculations.
                    """

                    # Note: May be type issues somewhere resulting in Any types and slowing down code.

                    # total_sum += coefficient * n_coefficient * k_even_coefficient * besselk(ArbReal(n + 1), ArbReal(abs(x) * a)) * abs(x)^(n + 1)
                end
            end
        end
    end

    # Return final value obtained from double summation.
    return total_sum
end

"""
ℑχ_0(Ω::Float64, α::Float64, v::Float64, w::Float64, N::Int)

    Calculate the imaginary part of χ(Ω) in a zero temperature approximation (equation (16) in Devreese's et al.) for a given frequency Ω. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. N is the upper limit of a sum that is analytically infinite, however most of the information is encapuslated by N < 50.
"""
function ℑχ_0(Ω, α, v, w, N = 10)

    R = (v^2 - w^2) / (w^2 * v)
    coefficient = 2 / 3 * α * (v / w)^3

    total_sum = 0.0
    for n = 0:Int(N)
        total_sum +=
            -coefficient * sqrt(π) / (gamma(-n - 1 / 2) * gamma(n + 1)) * (-2 * R)^n /
            BigCombinatorics.doublefactorial(2 * n + 1) *
            abs(Ω - 1 - n * v)^(n + 1 / 2) *
            exp(-R * abs(Ω - 1 - n * v)) *
            (1 + sign(Ω - 1 - n * v))
    end
    return total_sum
end

"""
----------------------------------------------------------------------
The Moblity of the Polaron.
----------------------------------------------------------------------
"""

"""
polaron_mobility(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64, N::Int)

    Calculate the moblity μ(Ω) of the polaron at finite temperatues (equation (1) in Hellwarth 1999) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. N is the upper limit of a sum in ℑχ(Ω).
"""
function polaron_mobility(Ω, β, α, v, w)
    1 / ℑχ(Ω, β, α, v, w)
end


"""
polaron_mobility_zero(Ω::Float64, α::Float64, v::Float64, w::Float64, N::Int)

    Calculate the moblity μ(Ω) of the polaron at zero-temperatures (equation (1) in Hellwarth 1999) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. N is the upper limit of a sum in ℑχ(Ω).
"""
function polaron_mobility_zero(Ω, α, v, w)
    1 / ℑχ_0(Ω, α, v, w)
end
