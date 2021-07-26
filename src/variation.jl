# variation.jl

"""
variation(α::Float64; v = 7.0, w = 6.0)

    Calculate v and w variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
    This version uses the original athermal action (Feynman 1955).
    Returns v, w.
"""
function variation(α; v = 0.0, w = 0.0)

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
        Fminbox(BFGS()),
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
            Fminbox(BFGS()),
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
function variation(α, β; v = 0.0, w = 0.0)

    # Intial guess for v and w.
    if v == 0.0 || w == 0.0 # Default values to start with. Generates a random float between 1.0 and 11.0
        initial = sort(rand(2), rev=true) .* 4.0 .+ 1.0
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
        Fminbox(BFGS()),
    )

    # Get v and w values that minimise the free energy.
    v, w = Optim.minimizer(solution)

    # If optimisation does not converge or if v ≤ w, pick new random starting guesses for v and w between 1.0 and 11.0. Repeat until the optimisation converges with v > w.
    while Optim.converged(solution) == false || v <= w
        initial = sort(rand(2), rev=true) .* 4.0 .+ 1.0
        solution = Optim.optimize(
            Optim.OnceDifferentiable(f, initial; autodiff = :forward),
            lower,
            upper,
            initial,
            Fminbox(BFGS()),
        )
        v, w = Optim.minimizer(solution)
    end

    # Return variational parameters that minimise the free energy.
    return v, w
end

"Multiple branch variation"

function multi_variation(T, ϵ_optic, m_eff, volume, freqs_and_ir_activity; N = 1) # N number of v and w params

    setprecision(BigFloat, 32) # Speed up. Stops potential overflows.

    # Number of phonon branches.
    M = length(freqs_and_ir_activity[:, 1])

    # Initialise MxN matrices for v and w parameters. M is number of phonon branches. N is number of variational parameters (v & w) per branch.
    v_params = Matrix{Float64}(undef, M, N)
    w_params = Matrix{Float64}(undef, M, N)

    # Intial guess for v and w.
    initial = sort(rand(2 * N)) .* 4.0 .+ 1.0 # initial guess around 4 and ≥ 1.

    # Limits of the optimisation.
    lower = repeat([0.1], 2 * N) 
    upper = repeat([200.0], 2 * N)

    for j in 1:M # sum over phonon branches

        # Osaka Free Energy function to minimise.
        f(x) = multi_free_energy([x[2 * n] for n in 1:Int(N)], [x[2 * n - 1] for n in 1:Int(N)], T, ϵ_optic, m_eff, volume, freqs_and_ir_activity, j)

        # Use Optim to optimise the free energy function w.r.t v and w.
        solution = Optim.optimize(
            Optim.OnceDifferentiable(f, initial; autodiff = :forward),
            lower,
            upper,
            initial,
            Fminbox(BFGS()),
            # Optim.Options(time_limit = 10.0),
        )

        # Get v and w params that minimised free energy.
        var_params = Optim.minimizer(solution)

        # If v ≤ w quit as negative mass.
        # if any(sort([var_params[2 * n] for n in 1:Int(N)]) .<= sort([var_params[2 * n - 1] for n in 1:Int(N)]))
        #     return "v_i <= w_i"
        # end

        # Intialise next guess of v and w to be their values for the current phonon branch. (quicker convergence)
        initial = sort(var_params)

        # Update matrices for v and w parameters.
        v_params[j, :] .= [var_params[2 * n] for n in 1:N]
        w_params[j, :] .= [var_params[2 * n - 1] for n in 1:N]

        # Show current v and w that minimise jth phonon branch.
        println(var_params)
    end

    # Return variational parameters that minimise the free energy.
    return v_params, w_params
end