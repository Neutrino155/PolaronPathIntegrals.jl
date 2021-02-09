# variation.jl

"""
variation(α::Float64; v = 7.0, w = 6.0)

    Calculate v and w variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
    This version uses the original athermal action (Feynman 1955).
    Returns v, w.
"""
function variation(α; v = 7.0, w = 6.0)

    initial = [v, w]
    lower = [0.0, 0.0]
    upper = [100.0, 100.0]

    f(x) = free_energy(x[1], x[2], α)

    od = Optim.OnceDifferentiable(f, initial; autodiff = :forward)
    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff = :forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
    )

    if Optim.converged(solution) == false
        print("\tWARNING: Failed to converge to v, w soln? : ", Optim.converged(solution))
    end

    v, w = Optim.minimizer(solution)
end

function variation(α, β; v = 1.534, w = 1.235)

    initial = [v, w]
    lower = [0.001, 0.001]
    upper = [Inf, v]

    f(x) = free_energy(x[1], x[2], α, β)

    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff = :forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
    )

    v, w = Optim.minimizer(solution)
end
