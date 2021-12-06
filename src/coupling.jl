# coupling.jl

"""
Implementation of Frohlich's Hamiltonian and coupling parameter α applied to the Polaron model.

See Frohlich 1952:
https://www.jstor.org/stable/99166

"""

"""
frohlich_α(ϵ_optic::Float64, ϵ_static::Float64, freq::Float64, m_eff::Float64)

    Calculates the Frohlich alpha parameter, for a given dielectric constant, frequency (f) of phonon in Hertz, and effective mass (in units of the bare electron mass).

    See Equation 2.37.

"""
function frohlich_α(ϵ_optic, ϵ_static, freq, m_eff)
    ω = 2π * freq * 1e12
    α = 0.5 / (4 * π * ϵ_0) *
       (1 / ϵ_optic - 1 / ϵ_static) *
       (eV^2 / (ħ * ω)) *
       sqrt(2 * m_e * m_eff * ω / ħ)
    return α
end

"""
----------------------------------------------------------------------
Multiple Branch Frohlich Alpha
----------------------------------------------------------------------

This section of the code is dedicated to determining the partial dielectric electron-phonon coupling parameter, decomposed into ionic dielectric contributions from each phonon mode of the matieral. 
"""

"""
ϵ_ionic_mode(phonon_mode_freq::Float64, ir_activity::Float64, volume::Float64)

    Calculate the ionic contribution to the dielectric function for a given phonon mode.
    phonon_mode_freq is the frequency of the mode in THz.

     - ir_activity is the infra-red activity of the mode in e^2 amu^-1.
     - volume is the volume of the unit cell of the material in m^3.
"""
function ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume) # single ionic mode

    # Angular phonon frequency for the phonon mode (rad Hz)
    ω_j = 2π * phonon_mode_freq * 1e12 

    # Dielectric contribution from a single ionic phonon mode
    ϵ_mode = eV^2 * ir_activity / (3 * volume * ω_j^2 * amu)

    # Normalise ionic dielectric contribution with 1 / (4π ϵ_0) (NB: the 4π has been pre-cancelled)
    return ϵ_mode / ϵ_0 
end

"""
ϵ_total(freqs_and_ir_activity::Matrix{Float64}, volume::Float64)

    Calculate the total ionic contribution to the dielectric function from all phonon modes.

     - freqs_and_ir_activity is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e^2 amu^-1) in the second column.
     - volume is the volume of the unit cell of the material in m^3.
"""
function ϵ_total(freqs_and_ir_activity, volume) # total ionic contribution to dielectric

    # Extract phonon frequencies (THz)
    phonon_freqs = freqs_and_ir_activity[:, 1] 

    # Extra infra-red activities (e^2 amu^-1)
    ir_activity = freqs_and_ir_activity[:, 2]

    # Sum over all ionic contribution from each phonon mode
    total_ionic = 0.0
    for (f, r) in zip(phonon_freqs, ir_activity)
        total_ionic += ϵ_ionic_mode(f, r, volume) 
    end

    return total_ionic
end

"""
effective_freqs(freqs_and_ir_activity::Matrix{Float64}, num_var_params::Integer)

    Generates a matrix of effective phonon modes with frequencies and infra-red activities derived from a larger matrix using the Principal Component Analysis (PCA) method.

     - freqs_and_ir_activity: is a matrix containing the phonon mode frequencies (in THz) in the first column and the infra-red activities (in e^2 amu^-1) in the second column.
     - num_var_params: is the number of effective modes required (which needs to be less than the number of modes in freqs_and_ir_activity)

    *** POSSIBLY REDUNDANT ***
"""
function effective_freqs(freqs_and_ir_activity, num_var_params) # PCA Algorithm

    # Check that the number of effective modes is less than the number of actual phonon modes.
    if num_var_params >= size(freqs_and_ir_activity)[1]

        println("The number of effective phonon modes has to be less than the total number of phonon modes.")

    else

        # Centralise data by subtracting the columnwise mean
        standardized_matrix = freqs_and_ir_activity' .- mean(freqs_and_ir_activity', dims = 2) 

        # Calculate the covariance matrix S' * S. Matrix size is (n - 1) x (n - 1) for number of params (here n = 2)
        covariance_matrix = standardized_matrix' * standardized_matrix 

        # Extract eigenvectors of the covariance matrix
        eigenvectors = eigvecs(covariance_matrix) 

        # Project the original data along the covariance matrix eigenvectors and undo the centralisation
        reduced_matrix = standardized_matrix[:, 1:num_var_params] * eigenvectors[1:num_var_params, 1:num_var_params] * 
        eigenvectors[1:num_var_params, 1:num_var_params]' .+ mean(freqs_and_ir_activity', dims = 2)

        # Resultant matrix is positive definite and transposed.
        return abs.(reduced_matrix')
    end
end

"""
frohlich_α_j(ϵ_optic::Float64, ϵ_ionic::Float64, ϵ_total::Float64, phonon_mode_freq::Float64, m_eff::Float64)

    Calculates the partial dielectric electron-phonon coupling parameter for a given longitudinal optical phonon mode. This decomposes the original Frohlich alpha coupling parameter (defined for a single phonon branch) into contributions from multiple phonon branches.

     - ϵ_optic is the optical dielectric constant of the material.
     - ϵ_ionic is the ionic dielectric contribution from the phonon mode.
     - ϵ_total is the total ionic dielectric contribution from all phonon modes of the material.
     - phonon_mode_freq is the frequency of the phonon mode (THz).
     - m_eff is the band mass of the electron (in units of electron mass m_e)
"""
function frohlich_α(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff) 

    # The Rydberg energy unit
    Ry = eV^4 * me / (2 * ħ^2)

    # Angular phonon frequency for the phonon mode (rad Hz).
    ω = 2π * 1e12 * phonon_mode_freq 

    # The static dielectric constant. Calculated here instead of inputted so that ionic modes are properly normalised.
    ϵ_static = ϵ_total + ϵ_optic

    # The contribution to the electron-phonon parameter from the currrent phonon mode. 1 / (4π ϵ_0) is the dielectric normalisation.
    α_j  = (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0) / (ϵ_optic * ϵ_static)

    return α_j
end