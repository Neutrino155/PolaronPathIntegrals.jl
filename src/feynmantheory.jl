#feynmantheory.jl

"""
Implementation of Feynman's original variational technique applied to the Polaron model.

See Feynman 1955:
http://dx.doi.org/10.1103/PhysRev.97.660
"""

# Equation 31: The <|X(t) - X(s)|^{-1}> * exp(-|t-w|) effective action.
A_integrand(v, w, τ) = (w^2 * τ + (v^2 - w^2) / v * (1 - exp(-v * τ)))^(-0.5) * exp(-τ)

A(v, w, α) = π^(-0.5) * α * v * QuadGK.quadgk(τ -> A_integrand(v, w, τ), 0, Inf)[1]

# Equation 33: Lowest Free energy E = -B - A where B = -3/(4v)*(v-w)^2.
feynman_free_energy(v, w, α) = (3 / (4 * v)) * (v - w)^2 - A(v, w, α)

"""
feynman_variation(α::Float64; v = 7.0, w = 6.0)

    Calculate v and w variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
    This version uses the original athermal action (Feynman 1955).
    Returns v, w.
"""
function feynman_variation(α; v = 7.0, w = 6.0)
    initial = [v, w]
    lower = [0.0, 0.0]
    upper = [100.0, 100.0]

    f(x) = E(x[1], x[2], α)
    od = Optim.OnceDifferentiable(f, initial; autodiff = :forward)
    solution = Optim.optimize(od, lower, upper, initial, Fminbox(BFGS()))
    v, w = Optim.minimizer(solution)

    return v, w
end
