using ArbNumerics
using SpecialFunctions
using QuadGK
using DSP
using ApproxFun
using BigCombinatorics
using Plots
plotly()
Plots.PlotlyBackend()
# r = range(7* pi, stop = 10 * pi, length = 1000)
# signal = [sinh(0.5 * BigFloat(x) - 5)*cos(3 *x + 1.5)/BigFloat(x) for x in r]
# # s = [sinh(0.5 * BigFloat(x) - 5) * sin(3 *BigFloat(x) + 1.5)/BigFloat(x) for x in r]
# s = [(1 - cosh(0.5 * BigFloat(x) - 5)) * sin(3 *BigFloat(x) + 1.5) / x for x in r]
#
# a = hilbert(signal)
# b = hilbert(s)
# p = plot(r, real(a))
# plot!(r, imag(b))
# display(p)

function ℜχ(Ω, β, α, v, w)

    # Set arguments to BigFloat precision. Without this the calculations break down due to large values of hyperbolic functions.

    Ω = BigFloat(Ω)
    β = BigFloat(β)
    α = BigFloat(α)
    v = BigFloat(v)
    w = BigFloat(w)

    # M(n, x) = 2 * (2 / x)^n * QuadGK.quadgk(t -> sin(x * t) * (1 + t^2)^(-n - 1/2), BigFloat(0), Inf, maxevals=10^3, order=10)[1] / (sqrt(π) * gamma(1/2 - n))

    M(n, x) = QuadGK.quadgk(t -> sin(x * t) * (1 + t^2)^(-n - 1/2), BigFloat(0), Inf, maxevals=10^3, order=10)[1]

    # Initialise constants.

    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    # Coefficient independent of summation variables n & k.

    coefficient =
        2 * α / (3 * sqrt(π)) * (v / w)^3 * β^(3 / 2) * sinh(β * Ω / 2) / sinh(β / 2)

    # Initialise total value of double summation as a BigFloat.

    total_sum = BigFloat(0.0)
    next_current = BigFloat(0.0)
    n = 0

    while next_current >= total_sum * 1e-3 # Infinite summation from the Binomial expansion of (x^2 + a^2 - b * cos(v * x))^(-3/2). |b * cos(v * x) / (x^2 + a^2)| < 1 for v > 0 and β > 0. Here just limit summation to when next term adds negiglible amount to total sum.

        current = BigFloat(0.0)

        # Coefficient that depends on n summation variable.

        n_coefficient = b^n * BigCombinatorics.doublefactorial(2 * n + 1) / 2^(2 * n + 1)

        # Summation over expansion of (cos(v * x))^n splits into sums over odd n and even n.

        # Finite sum over k for even n from (cos(v * x))^n expansion.

        for k = 0:floor((n - 1) / 2)

            # Coefficient that depends on k summation variable.

            k_coefficient_inverse = factorial(k) * factorial(n - k)

            # Arguments of resultant cosines obtained from cos^n(vx) * cos(x) * cos(Ωx).

            y = [
                (Ω + 1 + v * (n - 2 * k)),
                (Ω + 1 - v * (n - 2 * k)),
                -(Ω - 1 + v * (n - 2 * k)),
                -(Ω - 1 - v * (n - 2 * k))
            ]

            for (x, k_coeff) in zip(y, k_coefficient_inverse)  # Iterate over cosine arguments for brevity.

                """
                For low-ish values of β (~ < 40) and Ω (~ < 50). Integral over cos(yx) / (x^2 + a^2)^(n + 3/2) can be identified with BesselK functions.
                """

                spec = M(n + 1, abs(x * a))

                current +=
                    spec * coefficient *
                    n_coefficient / k_coeff
            end
        end
        total_sum += current
        next_current = current
        n += 1
    end

    # Return final value obtained from double summation.
    # @show(total_sum)
    return total_sum
end

# function M(n, x)
#     p = 128
#     setextrabits(0)
#     setworkingprecision(ArbReal, bits = p)
#     n = ArbReal(n)
#     x = ArbReal(x)
#     total_sum = ArbReal(0.0)
#     next_current = ArbReal(0.0)
#     k = ArbReal(0)
#
#     while abs(next_current) >= abs(total_sum) * 1e-10
#
#         current = sign(ArbReal(x, bits = p)) * ArbReal(x / 2, bits = p)^ArbReal(2 * k + n, bits = p) / (ArbNumerics.gamma(ArbReal(k + 1, bits = p)) * ArbNumerics.gamma(ArbReal(n + k + 1, bits = p))) - ArbNumerics.ArbReal(x / 2, bits = p)^ArbReal(2 * k - n + 1, bits = p) / (ArbNumerics.gamma(ArbReal(k + 3/2, bits = p)) * ArbNumerics.gamma(ArbReal(-n + k + 3/2, bits = p)))
#
#         next_current = ArbReal(current, bits = p)
#         total_sum += ArbReal(next_current, bits = p)
#
#         if floor(log10(abs(ball(total_sum)[1]))) == floor(log10(abs(ball(total_sum)[2]))) + 30
#             @show(typeof(total_sum), floor(log10(abs(ball(total_sum)[1]))), floor(log10(abs(ball(total_sum)[2]))))
#             p *= 2
#             setworkingprecision(ArbReal, bits = p)
#             total_sum = ArbReal(0.0, bits = p)
#             next_current = ArbReal(0.0, bits = p)
#             k = ArbReal(-1)
#         end
#         k += ArbReal(1)
#     end
#     return total_sum
# end

function extra(Ω, β, α, v, w)

    # Set arguments to BigFloat precision. Without this the calculations break down due to large values of hyperbolic functions.

    Ω = BigFloat(Ω)
    β = BigFloat(β)
    α = BigFloat(α)
    v = BigFloat(v)
    w = BigFloat(w)

    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient =
        2 * α * (v / w)^3 * β^(3 / 2) / (3 * sqrt(π) * sinh(β / 2))

    I1 = quadgk(x -> (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2), BigFloat(0), BigFloat(β / 2), maxevals=10^4, order=7)[1]
    #
    # I2 = quadgk(x -> sin(Ω * x) * cos(x) / (a^2 + x^2 - b * cos(v * x))^(3 / 2), BigFloat(0), Inf, maxevals=10^4, order=7)[1]
    # @show(coefficient * (sinh(Ω * β / 2) * I2 + I1))

    # I2 = 0.0
    # for n in 0:10
    #     I2 +=  - 2 * sqrt(π) * (-b)^n * quadgk(x -> (cos(x) * sin(Ω * x) * (cos(v * x))^n) / (x^2 + a^2)^(n + 3/2),  BigFloat(0), Inf, maxevals=10^4, order=7)[1] / (SpecialFunctions.gamma(n + 1) * SpecialFunctions.gamma(-1/2 - n))
    # end

    I2 = 0.0
    for n in 0:10
        n_coefficient = - 2 * sqrt(π) * (-b / 2)^n / (gamma(n + 1) * gamma(-1/2 - n))
        for k in -1:floor((n - 1)/2)
            I2 += n_coefficient * (2 * gamma(n + 1) * quadgk(x -> (sin(Ω * x) * cos(x) * cos((n - 2 * k) * v * x)) / (x^2 + a^2)^(n + 3/2), BigFloat(0), Inf, maxevals=10^4, order=7)[1] / (gamma(k + 1) * gamma(n - k + 1)) + ((1 + (-1)^n) / 2) * 2^n * gamma(n/2 + 1/2) * quadgk(x -> (sin(Ω * x) * cos(x)) / (x^2 + a^2)^(n + 3/2), BigFloat(0.0), Inf, maxevals=10^4, order=7)[1] / (sqrt(π) * gamma(n/2 + 1)))
        end
    end

    f = coefficient * (sinh(Ω * β / 2) * I2 + I1)
    @show(Ω, f)
    return f
end
#
# setprecision(BigFloat, 64)
# # t = ℜχ(30, 100, 7, 5.8, 1.6)
# r = range(0.01, stop = 10.0, length = 100)
# i = [extra(j, 3, 5.0, 4.0, 2.1) for j in r]
# p = plot(r, i)
# display(p)

function M(n, x)
    p = 128
    setextrabits(0)
    setworkingprecision(ArbReal, bits = p)
    n = ArbReal(n, bits = p)
    x = ArbReal(Float64(x), bits = p)
    total_sum = ArbReal(0.0, bits = p)
    next_current = ArbReal(0.0, bits = p)
    k = ArbReal(0, bits = p)

    while abs(next_current) >= abs(total_sum) * 1e-3

        current = sign(ArbReal(x, bits = p)) * ArbReal(x / 2, bits = p)^ArbReal(2 * k + n, bits = p) / (ArbNumerics.gamma(ArbReal(k + 1, bits = p)) * ArbNumerics.gamma(ArbReal(n + k + 1, bits = p))) - ArbNumerics.ArbReal(x / 2, bits = p)^ArbReal(2 * k - n + 1, bits = p) / (ArbNumerics.gamma(ArbReal(k + 3/2, bits = p)) * ArbNumerics.gamma(ArbReal(-n + k + 3/2, bits = p)))

        next_current = ArbReal(current, bits = p)
        total_sum += ArbReal(next_current, bits = p)

        if floor(log10(abs(ball(total_sum)[1]))) == floor(log10(abs(ball(total_sum)[2]))) + 30
            p *= 2
            setworkingprecision(ArbReal, bits = p)
            total_sum = ArbReal(0.0, bits = p)
            next_current = ArbReal(0.0, bits = p)
            k = ArbReal(-1)
        end
        k += ArbReal(1)
    end
    return total_sum
end
# n = 1
# x = 10
# a = 10.234234
# @time c = M(n + 1, a * abs(x)) * sign(x) * sqrt(pi) * gamma(-1/2 - n) * (abs(x) / (2 * a))^(n + 1) / 2
# @time b = quadgk(t -> sin(x * t) / (t^2 + a^2)^(n + 3/2), BigFloat(0.0), Inf, maxevals=10^7, order=7)[1]
# @show(c, b)
β = 3
v = 4.0
w = 2.1
Ω = 5
R = (v^2 - w^2) / (w^2 * v)
a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
b = R * β / sinh(β * v / 2)
I1 = Fun(x -> (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2), 0..(β/2))
# I2 = quadgk(x -> (1 - cosh(Ω * x)) * cosh(x - β / 2) / (a^2 - β^2 / 4 + x * (β - x) - b * cosh(v * (x - β / 2)))^(3 / 2), 0, β/2, maxevals = 10^4, order = 7)
# @show(I1)
@show(cumsum(I1(β/2))[1])
# @show(I2)

# for p in 100:10:250
#     setextrabits(0)
#     setworkingprecision(ArbFloat, bits = p)
#     bessely = ArbNumerics.bessely(2.0, parse(ArbReal, "57.9668555791947111700112064539202300037282162"))
#     @show((p, bessely))
# end

# a = ArbNumerics.besselk(2.0, ArbFloat(1e-400))
# @show(a, typeof(a), ArbNumerics.NaN)
# if a == NaN
#     println("True")
# end
