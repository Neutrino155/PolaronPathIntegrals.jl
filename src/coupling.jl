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
    ω = 2 * π * freq
    α = 0.5 / (4 * π * ϵ_0) *
       (1 / ϵ_optic - 1 / ϵ_static) *
       (eV^2 / (ħ * ω)) *
       sqrt(2 * m_e * m_eff * ω / ħ)
    return α
end

"Multiple branches frohlich α"

function ϵ_ionic_mode(phonon_mode_freq, ir_activity, volume) # single ionic mode
    ω_j = 2π * phonon_mode_freq * 1e12 # angular phonon freq in Hz
    ϵ_mode = eV^2 * ir_activity / (3 * volume * ω_j^2 * amu) # single ionic mode
    return ϵ_mode / ϵ_0 # normalise with 1 / (4π ϵ_0)
end

function ϵ_total(freqs_and_ir_activity, volume) # total ionic contribution to dielectric
    phonon_freqs = freqs_and_ir_activity[:, 1] 
    ir_activity = freqs_and_ir_activity[:, 2]
    result = 0.0
    for (f, r) in zip(phonon_freqs, ir_activity)
        result += ϵ_ionic_mode(f, r, volume) # sum over all ionic contributions
    end
    return result
end

function effective_freqs(freqs_and_ir_activity, num_var_params) #PCA Algorithm
    standardized_matrix = freqs_and_ir_activity' .- mean(freqs_and_ir_activity', dims = 2) # centralise data by subtracting columnwise mean
    covariance_matrix = standardized_matrix' * standardized_matrix # has 1 / (n - 1) for number of params n = 2
    eigenvectors = eigvecs(covariance_matrix) # eigenvectors to project data along
    reduced_matrix = # project data along eigenvectors and undo centralisation
    standardized_matrix[:, 1:num_var_params] * eigenvectors[1:num_var_params, 1:num_var_params] * 
    eigenvectors[1:num_var_params, 1:num_var_params]' .+ mean(freqs_and_ir_activity', dims = 2)
    return abs.(reduced_matrix')
end

function frohlich_α_j(ϵ_optic, ϵ_ionic, ϵ_total, phonon_mode_freq, m_eff) # Frohlich alpha decomposed into phonon branch contributions
    Ry = eV^4 * me / (2 * ħ^2) # Rydberg energy
    ω = 2π * 1e12 * phonon_mode_freq # angular phonon freq (Hz)
    ϵ_static = ϵ_total + ϵ_optic # static dielectric. Calculate here instead of input so that ionic modes properly normalised.
    return (m_eff * Ry / (ħ * ω))^(1 / 2) * ϵ_ionic / (4π * ϵ_0) / (ϵ_optic * ϵ_static) # 1 / (4π ϵ_0) dielectric normalisation
end
