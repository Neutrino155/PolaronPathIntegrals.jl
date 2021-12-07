# optical_absorption.jl

"""
----------------------------------------------------------------------
Polaron absorption coefficient Γ(Ω).
----------------------------------------------------------------------
"""

"""
optical_absorption(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64; rtol = 1e-3)

    Calculate the absorption coefficient Γ(Ω) for the polaron at at finite temperatures (equation (11a) in [1]) for a given frequency Ω.  β is thermodynamic beta. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral in the memory function.  

    [1] Devreese, J., De Sitter, J., & Goovaerts, M. (1972). Optical Absorption of Polarons in the Feynman-Hellwarth-Iddings-Platzman Approximation. Physical Review B, 5(6), 2367–2381. doi:10.1103/physrevb.5.2367 

"""
function optical_absorption(Ω, β, α, v, w; rotl = 1e-3)
    real(complex_conductivity(Ω, β, α, v, w; rtol = rtol))
end

"""
----------------------------------------------------------------------
The Complex Impedence and Conductivity of the Polaron.
----------------------------------------------------------------------
"""

"""
polaron_complex_impedence(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the complex impedence Z(Ω) of the polaron at finite temperatures for a given frequency Ω (equation (41) in FHIP 1962 [1]). β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral in the memory function. 
"""
function polaron_complex_impedence(Ω, β, α, v, w; ω = 0.0, rtol = 1e-3)
	return -im * 2π * Ω / (sum(ω) * length(ω)) + im * polaron_memory_function(Ω, β, α, v, w; ω = ω, rtol = rtol)
end

"""
polaron_complex_conductivity(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the complex conductivity σ(Ω) of the polaron at finite temperatures for a given frequency Ω (equal to 1 / Z(Ω) with Z the complex impedence). β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error tolerance for the QuadGK integral in the memory function. 
"""
function polaron_complex_conductivity(Ω, β, α, v, w; ω = 0.0, rtol = 1e-3)
	return 1 / polaron_complex_impedence(Ω, β, α, v, w; ω = ω, rtol = rtol)
end

"""
function multi_conductivity(Ω::Float64, β::Array{Float64}(undef, 1), α::Array{Float64}(undef, 1), v::Array{Float64}(undef, 1), w::Array{Float64}(undef, 1), ω::Array{Float64}(undef, 1), m_eff::Float64)

    Calculate polaron complex conductivity inclusive of multiple phonon branches j, each with angular frequency ω[j] (rad THz).

     - Ω is the frequency (THz) of applied electric field.
     - β is an array of reduced thermodynamic betas, one for each phonon frequency ω[j]. 
     - α is an array of decomposed Frohlich alphas, one for each phonon frequency ω[j]. 
     - v is an one-dimensional array of the v variational parameters.
     - w is an one-dimensional array of the w variational parameters.
     - m_eff is the is the conduction band mass of the particle (typically electron / hole, in units of electron mass m_e).
"""
function multi_conductivity(Ω, β, α, v, w, ω, m_eff)

	# FHIP1962, page 1009, eqn (36).
    S(t, β) = cos(t - 1im * β / 2) / sinh(β / 2) / D_j(-1im * t, β, v, w)^(3 / 2)

    # FHIP1962, page 1009, eqn (35a).
    integrand(t, β, Ω) = (1 - exp(1im * 2π * Ω * t)) * imag(S(t, β))

    impedence = 0.0

    # Sum over the phonon modes.
	for j in 1:length(ω)

        # print out the current photon frequency and phonon mode frequency (THz).
		# println("Photon frequency = $Ω, Phonon mode frequency = $(ω[j] / 2π)")

        # Add the contribution to the complex impedence from the `jth` phonon mode.
		impedence += -1im * 2π * Ω / length(ω) + 1im * 2 * α[j] * ω[j]^2 * quadgk(t -> integrand(t, β[j], Ω / ω[j]), 0.0, Inf)[1] / (3 * √π * 2π * Ω)
	end

    # Conductivity is the reciprocal of the impedence. 
	conductivity = 1 / impedence * eV * 100^2 / (m_eff * me * 1e12)

    # Print out the value of the complex conductivity.
    # println("Multiple complex conductivity: ", conductivity, " cm^2/Vs")

    return conductivity
end