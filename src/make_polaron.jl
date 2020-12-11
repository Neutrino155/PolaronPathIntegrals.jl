# make_polaron.jl

function make_polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff; temp = 300.0, efield_freq = 0.001)
    ω = 2 * π * phonon_freq
    Ω = efield_freq
    T = temp
    α = frohlich_α(ϵ_optic, ϵ_static, phonon_freq, m_eff)
    if typeof(temp) == Float64
        if temp == 0.0
            β = :∞
            v, w = feynman_variation(α)
            F = feynman_free_energy(v, w, α)
            if typeof(Ω) != Float64
                μ = []
                Γ = []
                for i in Ω
                    append!(μ, 100^2 * eV * polaron_mobility_zero(i, α, v, w) / (ω * m_eff * m_e))
                    append!(Γ, optical_absorption_zero(i, α, v, w, sqrt(ϵ_optic)))
                end
            else
                μ = 100^2 * eV * polaron_mobility_zero(Ω, α, v, w) / (ω * m_eff * m_e)
                Γ = optical_absorption_zero(Ω, α, v, w, sqrt(ϵ_optic))
            end
        elseif temp != 0.0
            β = ħ * ω / (k_B * temp)
            v, w = singlemode_variation(α, β)
            F = osaka_free_energy(v, w, β, α)
            if typeof(Ω) != Float64
                μ = []
                Γ = []
                for i in Ω
                    append!(μ, 100^2 * eV * polaron_mobility(i, β, α, v, w) / (ω * m_eff * m_e))
                    append!(Γ, optical_absorption(i, β, α, v, w))
                end
            else
                μ = 100^2 * eV * polaron_mobility(Ω, β, α, v, w) / (ω * m_eff * m_e)
                Γ = optical_absorption(Ω, β, α, v, w)
            end
        end
    elseif typeof(temp) != Float64
        β = [ħ * ω / (k_B * i) for i in temp]
        v = [singlemode_variation(α, i)[1] for i in β]
        w = [singlemode_variation(α, i)[2] for i in β]
        F = [osaka_free_energy(i, j, k, α) for (i, j, k) in zip(v, w, β)]
        μ = []
        Γ = []
        for (i, j, k) in zip(v, w, β)
            append!(μ, 100^2 * eV * polaron_mobility(Ω, k, α, i, j) / (ω * m_eff * m_e))
            append!(Γ, optical_absorption(Ω, k, α, i, j))
        end
    end
    return Polaron(α, T, β, v, w, F, Ω, μ, Γ)
end
