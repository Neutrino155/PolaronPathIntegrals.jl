# variation.jl

"""
variation(α::Float64; v = 7.0, w = 6.0)

    Calculate v and w variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
    This version uses the original athermal action (Feynman 1955).
    Returns v, w.
"""
function variation(α; v = 5.0, w = 4.0)

    # Intial guess for v and w.
    if v == 0.0 || w == 0.0 # Default values to start with. Generates a random float between 1.0 and 11.0
        initial = sort(rand(2), rev=true) .* 10.0 .+ 1.0
    else
        initial = [v, w]
    end

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [Inf, Inf]

    # Osaka Free Energy function to minimise.
    f(x) = free_energy(x[1], x[2], α)

    # Use Optim to optimise the free energy function w.r.t v and w.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff = :forward),
        lower,
        upper,
        initial,
        Fminbox(LBFGS()),
    )

    # Get v and w values that minimise the free energy.
    v, w = Optim.minimizer(solution)

    # If optimisation does not converge or if v ≤ w, pick new random starting guesses for v and w between 1.0 and 11.0. Repeat until the optimisation converges with v > w.
    while Optim.converged(solution) == false || v <= w
        initial = sort(rand(2), rev=true) .* 10.0 .+ 1.0
        solution = Optim.optimize(
            Optim.OnceDifferentiable(f, initial; autodiff = :forward),
            lower,
            upper,
            initial,
            Fminbox(LBFGS()),
        )
        v, w = Optim.minimizer(solution)
    end

    # Return variational parameters that minimise the free energy.
    return v, w
end

"""
variation(α::Float64, β::Float64; v = 3.0, w = 2.0)

    Calculate v and w variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling and β reduced thermodynamic beta.
    This version uses temperature dependent action (Osaka 1959).
    Returns v, w.
"""
function variation(α, β, v = 3.0, w = 2.0)

    # Intial guess for v and w.
    if v == 0.0 || w == 0.0 # Default values to start with. Generates a random float between 1.0 and 11.0
        initial = sort(rand(2), rev=true) .* 10.0 .+ 1.0
    else
        initial = [v, w]
    end

    # Limits of the optimisation.
    lower = [0.0, 0.0]
    upper = [Inf, Inf]

    # Osaka Free Energy function to minimise.
    f(x) = free_energy(x[1], x[2], α, β)

    # Use Optim to optimise the free energy function w.r.t v and w.
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff = :forward),
        lower,
        upper,
        initial,
        Fminbox(LBFGS()),
    )

    # Get v and w values that minimise the free energy.
    v, w = Optim.minimizer(solution)

    # If optimisation does not converge or if v ≤ w, pick new random starting guesses for v and w between 1.0 and 11.0. Repeat until the optimisation converges with v > w.
    while Optim.converged(solution) == false || v <= w
        initial = sort(rand(2), rev=true) .* 10.0 .+ 1.0
        solution = Optim.optimize(
            Optim.OnceDifferentiable(f, initial; autodiff = :forward),
            lower,
            upper,
            initial,
            Fminbox(LBFGS()),
        )
        v, w = Optim.minimizer(solution)
    end

    # Return variational parameters that minimise the free energy.
    return v, w
end
