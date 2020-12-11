#hellwarththeory.jl

"""
Implementation of Hellwarth's extended treatment of Osaka's finite temperature version of Feynman's variational technique applied to the Polaron model; generalising it from the case of a single phonon mode to the case of multiple phonon modes.

And Hellwarth 1999:
https://doi.org/10.1103/PhysRevB.60.299

"""

function singlemode_variation(α, β; v = 7.2, w = 6.5)

    initial = [v, w]
    lower = [0.001, 0.001]
    upper = [Inf, v]

    f(x) = osaka_free_energy(x[1], x[2], β, α)

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
    return v, w
end
