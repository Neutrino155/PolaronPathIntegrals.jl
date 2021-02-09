# make_polaron.jl

# Example: MAPI: make_polaron(4.5, 24.1, 2.25e12, 0.12)

function make_polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff; temp = 300.0, efield_freq = 0.0)

    # Collect data.
    ω = 2 * π * phonon_freq
    Ω = efield_freq
    T = temp
    α = frohlich_α(ϵ_optic, ϵ_static, phonon_freq, m_eff)

    # Prepare empty arrays.
    β = []
    v = []
    w = []
    F = []
    μ = []
    Γ = []

    for t in T
        if t == 0.0
            append!(β, Inf)
            V, W = variation(α)
            append!(v, V)
            append!(w, W)
            append!(F, free_energy(V, W, α))

            μ_T = []
            Γ_T = []
            for f in Ω
                f = abs(f)
                append!(μ_T, 100^2 * eV * polaron_mobility(f, α, V, W, 1e4) / (ω * m_eff * m_e))
                append!(Γ_T, optical_absorption(f, α, V, W, 1e4)) / (c * ϵ_0 * sqrt(ϵ_optic))
            end
            append!(μ, [μ_T])
            append!(Γ, [Γ_T])

        elseif t > 0.0
            b = ħ * ω / (k_B * t)
            V, W = variation(α, b)
            @show(V, W, t)
            append!(β, b)
            append!(v, V)
            append!(w, W)
            @show(free_energy(V, W, α, b))
            append!(F, free_energy(V, W, α, b))

            μ_T = []
            Γ_T = []
            for f in Ω
                f = abs(f)
                append!(μ_T, 100^2 * eV * polaron_mobility(f, α, V, W, b) / (ω * m_eff * m_e))
                append!(Γ_T, optical_absorption(f, α, V, W, b))
            end
            append!(μ, [μ_T])
            append!(Γ, [Γ_T])

        else
            println("Temperature must be either zero or positive.")
        end
    end

    return Polaron(α, T, β, v, w, F, Ω, μ, Γ)
end
