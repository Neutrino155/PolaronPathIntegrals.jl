function polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff; temp = 0.0, field_freq = 0.0)
    ω = 2 * π * phonon_freq
    Ω = field_freq
    α = frohlich_α(ϵ_optic, ϵ_static, phonon_freq, m_eff)
    if temp == 0.0
        v, w = feynman_variation(α; v = 7.2, w = 6.5)
        F = feynman_free_energy(v, w, α)
        μ = ω * m_eff * polaron_mobility_zero(Ω, α, v, w) / eV
    else temp != 0.0
        β = ħ * ω / (k_B * temp)
        v, w = singlemode_variation(α, β; v = 7.2, w = 6.5)
        F = osaka_free_energy(v, w, β, α)
        μ = ω * m_eff * polaron_mobility(Ω, β, α, v, w) / eV
    end
    return α, β, v, w, F, μ
end
