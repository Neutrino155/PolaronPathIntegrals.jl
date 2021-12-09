# free_energy.jl

"""
Implementation of Feynman's original variational technique applied to the Polaron model.

See Feynman 1955:
http://dx.doi.org/10.1103/PhysRev.97.660
"""

# Equation 31: The <|X(t) - X(s)|^{-1}> * exp(-|t-w|) effective action.
A_integrand(x, v, w) = (abs(w^2 * x + (v^2 - w^2) / v * (1 - exp(-v * x))))^(-0.5) * exp(-x)

A(v, w, α) = π^(-0.5) * α * v * QuadGK.quadgk(x -> A_integrand(x, v, w), 0, Inf)[1]

# Equation 33: Lowest Free energy E = -B - A where B = -3/(4v)*(v-w)^2.
free_energy(v, w, α; ω = 1.0) = ((3 / (4 * v)) * (v - w)^2 - A(v, w, α)) * ω

"""
    Implementation of Osaka's extended treatment of Feynman's variational technique applied to the Polaron model; generalising it from the case at O^{∘}K to the case at finite temperature.

    See Osaka 1959:
    Progress of Theoretical Physics, Vol. 22, No.3, September 1959

    And Hellwarth 1999:
    https://doi.org/10.1103/PhysRevB.60.299

"""

# Equation 62d in Hellwarth.
Y(x, v, β) = 1 / (1 - exp(-v * β)) * (1 + exp(-v * β) - exp(-v * x) - exp(v * (x - β)))

# Integrand of Equation 62c in Hellwarth.
A_integrand(x, v, w, β) =  (exp(β - x) + exp(x)) / sqrt(abs(w^2 * x * (1 - x / β) + Y(x, v, β) * (v^2 - w^2) / v))

# Equation 62c in Hellwarth.
A(v, w, α, β) = α * v / (sqrt(π) * (exp(BigFloat(β)) - 1)) * QuadGK.quadgk(x -> A_integrand(x, v, w, β), BigFloat(0.0), BigFloat(β / 2))[1]

# Equation 62b in Hellwarth. Equation 20 in Osaka.
B(v, w, β) = 3 / β * (log(v / w) - 1 / 2 * log(2 * π * BigFloat(β)) - log(sinh(v * BigFloat(β) / 2) / sinh(w * BigFloat(β) / 2)))

# Equation 62e in Hellwarth. Equation 17 in Osaka.
C(v, w, β) = 3 / 4 * (v^2 - w^2) / v * (coth(v * BigFloat(β) / 2) - 2 / (v * β))

# Equation 62a in Hellwarth. In paragraph below Equation 22 in Osaka; has extra 1/β due to different definition of A, B & C.
function free_energy(v, w, α, β; ω = 1.0)
    setprecision(BigFloat, 64)
    a = A(v, w, α, β[1])
    b = B(v, w, β[1])
    c = C(v, w, β[1])
    -(a + b + c) * ω
end

"""
----------------------------------------------------------------------
Multiple Branch Polaron Free Energy
----------------------------------------------------------------------

This section of the code is dedicated to calculating the polaron free energy, generalised from Osaka's expression to the case where multiple phonon modes are present in the material.
"""

"""
κ_i(i::Integer, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Calculates the spring-constant coupling the electron to the `ith' fictitious mass that approximates the exact electron-phonon interaction with a harmonic coupling to a massive fictitious particle. NB: Not to be confused with the number of physical phonon branches; many phonon branches could be approximated with one or two etc. fictitious masses for example. The number of fictitious mass does not necessarily need to match the number of phonon branches.

     - i enumerates the current fictitious mass.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function κ_i(i, v, w)
    κ = v[i]^2 - w[i]^2
    if length(v) > 1
        for j in 1:length(v)
            if j != i
                κ *= (v[j]^2 - w[i]^2) / (w[j]^2 - w[i]^2)
            end
        end
    end
    return κ
end

"""
h_i(i::Integer, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Calculates the normal-mode (the eigenmodes) frequency of the coupling between the electron and the `ith' fictitious mass that approximates the exact electron-phonon interaction with a harmonic coupling to a massive fictitious particle. NB: Not to be confused with the number of physical phonon branches; many phonon branches could be approximated with one or two etc. fictitious masses for example. The number of fictitious mass does not necessarily need to match the number of phonon branches.

     - i enumerates the current fictitious mass.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function h_i(i, v, w)
    h = v[i]^2 - w[i]^2
    if length(v) > 1
        for j in 1:length(v)
            if j != i
                h *= (w[j]^2 - v[i]^2) / (v[j]^2 - v[i]^2)
            end
        end
    end
    return h
end

"""
C_ij(i::Integer, j::Integer, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Calculates the element to the coupling matrix C_ij (a generalisation of Feynman's `C` coupling variational parameter) between the electron and the `ith' and `jth' fictitious masses that approximates the exact electron-phonon interaction with a harmonic coupling to a massive fictitious particle. NB: Not to be confused with the number of physical phonon branches; many phonon branches could be approximated with one or two etc. fictitious masses for example. The number of fictitious mass does not necessarily need to match the number of phonon branches.

     - i, j enumerate the current fictitious masses under focus (also the index of the element in the coupling matrix C)
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function C_ij(i, j, v, w)
    C = w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
    return C
end

"""
D_j(τ::Float64, β::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Calculates the recoil function that approximates the exact influence (recoil effects) of the phonon bath on the electron with the influence of the harmonicly coupled fictitious masses on the electron. It appears in the exponent of the intermediate scattering function.

     - τ is the imaginary time variable.
     - β is the reduced thermodynamic temperature ħ ω_j / (kB T) associated with the `jth` phonon mode.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function D_j(τ, β, v, w)
    D = τ * (1 - τ / β)
    for i in 1:length(v)
        if v[i] != w[i]
        D += (h_i(i, v, w) / v[i]^2) * (2 * sinh(v[i] * τ / 2) * sinh(v[i] * (β - τ) / 2) / (v[i] * sinh(v[i] * β / 2)) - τ * (1 - τ / β))
        end
    end
    return D
end

"""
D_j(τ::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Calculates the recoil function that approximates the exact influence (recoil effects) of the phonon bath on the electron with the influence of the harmonicly coupled fictitious masses on the electron. It appears in the exponent of the intermediate scattering function. This function works at zero temperature.

     - τ is the imaginary time variable.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function D_j(τ, v, w)
    D = τ
    for i in 1:length(v)
        D += (h_i(i, v, w) / v[i]^2) * ((1 - exp(-v[i] * τ)) / v[i] - τ)
    end
    return D
end

"""
B_j(α::Float64, β::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Generalisation of the B function from Equation 62c in Biaggio and Hellwarth []. This is the expected value of the exact action <S_j> taken w.r.t trial action, given for the `jth` phonon mode.

     - α is the partial dielectric electron-phonon coupling parameter for the `jth` phonon mode.
     - β is the reduced thermodynamic temperature ħ ω_j / (kB T) associated with the `jth` phonon mode.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function B_j(α, β, v, w)
    B_integrand(τ) = cosh(β / 2 - abs(τ)) / (sinh(β / 2) * sqrt(abs(D_j(abs(τ), β, v, w))))
    B = α / √π * quadgk(τ -> B_integrand(τ), 0.0, β / 2)[1]
    return B
end

"""
B_j(α::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1))

    Generalisation of the B function from Equation 62c in Biaggio and Hellwarth [] to zero temperature. This is the expected value of the exact action <S_j> taken w.r.t trial action, given for the `jth` phonon mode.

     - α is the partial dielectric electron-phonon coupling parameter for the `jth` phonon mode.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
"""
function B_j(α, v, w)
    B_integrand(τ) = exp(-abs(τ)) / sqrt(abs(D_j(abs(τ), v, w)))
    B = α / √π * quadgk(τ -> B_integrand(τ), 0.0, Inf)[1]
    return B
end

"""
C_j(β::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), n::Float64)

    Generalisation of the C function from Equation 62e in Biaggio and Hellwarth []. This is the expected value of the trial action <S_0> taken w.r.t trial action.

     - β is the reduced thermodynamic temperature ħ ω_j / (kB T) associated with the `jth` phonon mode.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - n is the total number of phonon modes.
"""
function C_j(β, v, w, n)

    s = 0.0

    # Sum over the contributions from each fictitious mass.
    for i in 1:length(v)
        for j in 1:length(v)
            s += C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2)  - 2 / (β * v[j]))
        end
    end

    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    return 3 * s / n
end

"""
C_j(v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), n::Float64)

    Generalisation of the C function from Equation 62e in Biaggio and Hellwarth [] but to zero temperaure. This is the expected value of the trial action <S_0> taken w.r.t trial action.

     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - n is the total number of phonon modes.
"""
function C_j(v, w, n)

    s = 0.0

    # Sum over the contributions from each fictitious mass.
    for i in 1:length(v)
        for j in 1:length(v)
            s += C_ij(i, j, v, w) / (v[j] * w[i])
        end
    end

    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    return 3 * s / n
end

"""
A_j(β::Float64, v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), n::Float64)

    Generalisation of the A function from Equation 62b in Biaggio and Hellwarth []. This is the Helmholtz free energy of the trial model.

     - β is the reduced thermodynamic temperature ħ ω_j / (kB T) associated with the `jth` phonon mode.
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - n is the total number of phonon modes.
"""
function A_j(β, v, w, n)

    s = -log(2π * β) / 2

    # Sum over the contributions from each fictitious mass.
    for i in 1:length(v)
        if v[i] != w[i]
            s += log(v[i] / w[i]) - log(sinh(v[i] * β / 2) / sinh(w[i] * β / 2))
        end
    end

    # Divide by the number of phonon modes to give an average contribution per phonon mode.
    3 / β * s / n
end

"""
A_j(v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), n::Float64)

    Generalisation of the A function from Equation 62b in Biaggio and Hellwarth [] but to zero temperature. This is the ground-state energy of the trial model.

     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - n is the total number of phonon modes.
"""
function A_j(v, w, n)
    s = 0.0
    for i in 1:length(v)
        s += v[i] - w[i]
    end
    return -3 * s / (2 * n)
end

"""
free_energy(v_params::Array{Float64}(undef, 1), w_params::Array{Float64}(undef, 1), T::Float64, ϵ_optic::Float64, m_eff::Float64, volume::Float64, freqs_and_ir_activity::Matrix{Float64})

    Calculates the Helmholtz free energy of the polaron for a material with multiple phonon branches. This generalises Osaka's free energy expression (below Equation (22) in []).

     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - T is the environment temperature in kelvin (K).
     - ϵ_optic is the optical dielectric constant of the material.
     - m_eff is the band mass of the electron (in units of electron mass m_e) .
     - volume is the volume of the unit cell of the material in m^3.
     - freqs_and_ir_activity is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e^2 amu^-1) in the second column.
"""
function free_energy(v, w, α::Array, β::Array; ω = 1.0)

    # Total number of phonon modes / branches.
    num_of_modes = length(ω)

    F = 0.0

    # Sum over the phonon modes.
	for j in 1:num_of_modes

        # Add contribution to the total free energy from the phonon mode.
		F += -(B_j(α[j], β[j], v, w) + C_j(β[j], v, w, num_of_modes) + A_j(β[j], v, w, num_of_modes)) * ω[j]
        
        # Prints out the frequency, reduced thermodynamic temperature, ionic dielectric and partial coupling for the phonon mode.
        # println("Free energy: Phonon freq = ", phonon_freqs[j], " | β = ", β_j, " | ϵ_ionic = ", ϵ_ionic_j, " | α_j = ", α_j)
    end
	
    # print out the total polaron free energy from all phonon modes.
    # println("Total free energy: ", F * ħ / eV * 1e3, " meV")

    # Free energy in units of meV
    return F
end

"""
free_energy(v_params::Array{Float64}(undef, 1), w_params::Array{Float64}(undef, 1))

    Calculates the Helmholtz free energy of the polaron for a material with multiple phonon branches. This generalises Osaka's free energy expression (below Equation (22) in []).

     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - T is the environment temperature in kelvin (K).
     - ϵ_optic is the optical dielectric constant of the material.
     - m_eff is the band mass of the electron (in units of electron mass m_e) .
     - volume is the volume of the unit cell of the material in m^3.
     - freqs_and_ir_activity is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e^2 amu^-1) in the second column.
"""
function free_energy(v, w, α::Array; ω = 1.0)

    # Speed up. Stops potential overflows.
    setprecision(BigFloat, 32) 

    # Total number of phonon modes / branches.
    num_of_branches = length(ω)

    F = 0.0

    # Sum over the phonon modes.
	for j in 1:num_of_branches

        # Add contribution to the total free energy from the phonon mode.
		F += -(B_j(α[j], v, w) + C_j(v, w, num_of_branches) + A_j(v, w, num_of_branches)) * ω[j]
        
        # Prints out the frequency, reduced thermodynamic temperature, ionic dielectric and partial coupling for the phonon mode.
        # println("Free energy: Phonon freq = ", phonon_freqs[j], " | β = ", β_j, " | ϵ_ionic = ", ϵ_ionic_j, " | α_j = ", α_j)
    end
	
    # print out the total polaron free energy from all phonon modes.
    # println("Total free energy: ", F * ħ / eV * 1e3, " meV")

    # Free energy in units of meV
    return F
end