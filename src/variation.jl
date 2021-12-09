# variation.jl

"""
variation(α::Float64; v = 7.0, w = 6.0)

    Calculate v and w variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
    This version uses the original athermal action (Feynman 1955).
    Returns v, w.
"""
function variation(α; v = 0.0, w = 0.0, ω = 1.0)

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
    f(x) = free_energy(x[1], x[2], α; ω = ω)

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

    # Return variational parameters that minimise the free energy.
    return v, w
end

"""
variation(α::Float64, β::Float64; v = 3.0, w = 2.0)

    Calculate v and w variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling and β reduced thermodynamic beta.
    This version uses temperature dependent action (Osaka 1959).
    Returns v, w.
"""
function variation(α, β; v = 0.0, w = 0.0, ω = 1.0)

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
    f(x) = free_energy(x[1], x[2], α, β; ω = ω)

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

    # Return variational parameters that minimise the free energy.
    return v, w
end

"""
multi_variation(T::Float64, ϵ_optic::Float64, m_eff::Float64, volume::Float64, freqs_and_ir_activity::Matrix{Float64}; initial_vw::{Bool, Array{Float64}(undef, 1)}, N::Integer)

    Minimises the multiple phonon mode free energy function for a set of v_p and w_p variational parameters.
    The variational parameters follow the inequality: v_1 > w_1 > v_2 > w_2 > ... > v_N > w_N.

     - T is the environment temperature in kelvin (K).
     - ϵ_optic is the optical dielectric constant of the material.
     - m_eff is the band mass of the electron (in units of electron mass m_e) .
     - volume is the volume of the unit cell of the material in m^3.
     - freqs_and_ir_activity is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e^2 amu^-1) in the second column.
     - initial_vw determines if the function should start with a random initial set of variational parameters (Bool input) or a given set of variational parameter values (one dimensional array).
     - N specifies the number of variational parameter pairs, v_p and w_p, to use in minimising the free energy.
"""
function variation(α::Array, β::Array; v = 0.0, w = 0.0, ω = 1.0, N = 1, T = 20.0) # N number of v and w params

    if N != length(v) != length(w)
        return error("The number of variational parameters v & w must be equal to N.")
    end

    # Use a random set of N initial v and w values.
    if v == 0.0 || w == 0.0
		# Intial guess for v and w parameters.
    	initial = sort(rand(2 * N), rev=true) .* 4.0 .+ 1.0 # initial guess around 4 and ≥ 1.
	else
        initial = vcat(v, w)
    end

    # Limits of the optimisation.
    lower = fill(0.0, 2 * N)
    upper = fill(100.0, 2 * N)

    # Print out the initial v and w values.
	# println("Initial guess: ", initial)

	# The multiple phonon mode free energy function to minimise.
	f(x) = free_energy([x[2 * n - 1] for n in 1:N], [x[2 * n] for n in 1:N], α, β; ω = ω)

	# Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
	solution = Optim.optimize(
		f,
		lower,
		upper,
		initial,
		SAMIN(),
		Optim.Options(f_reltol = 1e-3, x_reltol = 1e-3, iterations=10^6, show_trace = true, show_every = 50), # Set time limit for asymptotic convergence if needed.
	)

	# Extract the v and w parameters that minimised the free energy.
	var_params = sort(Optim.minimizer(solution), rev = true)

	# Separate the v and w parameters into one-dimensional arrays (vectors).
	v = [var_params[2 * n - 1] for n in 1:N]
	w = [var_params[2 * n] for n in 1:N]

	# Print the variational parameters that minimised the free energy.
	# println("Variational parameters: ", var_params)

    # Return the variational parameters that minimised the free energy.
    return v, w
end

function variation(α::Array; v = 0.0, w = 0.0, ω = 1.0, N = 1, T = 20) # N number of v and w params
 
    if N != length(v) != length(w)
        return error("The number of variational parameters v & w must be equal to N.")
    end

    # Use a random set of N initial v and w values.
    if v == 0.0 || w == 0.0
		# Intial guess for v and w parameters.
    	initial = sort(rand(2 * N), rev = true) .* 4.0 .+ 1.0 # initial guess around 4 and ≥ 1.
	else
        initial = vcat(v, w)
    end

    # Limits of the optimisation.
    lower = fill(0.0, 2 * N)
    upper = fill(Inf, 2 * N)

    # Print out the initial v and w values.
	# println("Initial guess: ", initial)

	# The multiple phonon mode free energy function to minimise.
	f(x) = free_energy([x[2 * n - 1] for n in 1:N], [x[2 * n] for n in 1:N], α; ω = ω)

	# Use Optim to optimise the free energy function w.r.t the set of v and w parameters.
	solution = Optim.optimize(
		Optim.OnceDifferentiable(f, initial; autodiff = :forward),
		lower,
		upper,
		initial,
		Fminbox(LBFGS()),
		Optim.Options(f_reltol = 1e-3, x_reltol = 1e-3, show_trace = true, show_every = 50), # Set time limit for asymptotic convergence if needed.
	)

	# Extract the v and w parameters that minimised the free energy.
	var_params = sort(Optim.minimizer(solution), rev = true)

	# Separate the v and w parameters into one-dimensional arrays (vectors).
	v = [var_params[2 * n - 1] for n in 1:N]
	w = [var_params[2 * n] for n in 1:N]

	# Print the variational parameters that minimised the free energy.
	# println("Variational parameters: ", var_params)

    # Return the variational parameters that minimised the free energy.
    return v, w
end