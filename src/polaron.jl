# polaron.jl

function make_polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff; temp = 0.0, field_freq = 0.001)
    ω = 2 * π * phonon_freq
    Ω = field_freq
    α = frohlich_α(ϵ_optic, ϵ_static, phonon_freq, m_eff)
    if temp == 0.0
        β = "Infinity"
        v, w = feynman_variation(α)
        F = feynman_free_energy(v, w, α)
        μ = ω * m_eff * polaron_mobility_zero(Ω, α, v, w) / eV
        Γ = optical_absorption_zero(Ω, α, v, w, sqrt(ϵ_optic))
    else temp != 0.0
        β = ħ * ω / (k_B * temp)
        v, w = singlemode_variation(α, β)
        F = osaka_free_energy(v, w, β, α)
        μ = eV * polaron_mobility(Ω, β, α, v, w) / (ω * m_eff * m_e)
        Γ = optical_absorption(Ω, β, α, v, w)
    end
    return Polaron(α, β, v, w, F, μ, Γ)
end
