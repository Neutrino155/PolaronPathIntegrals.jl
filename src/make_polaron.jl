# make_polaron.jl

# Example: MAPI: make_polaron(4.5, 24.1, 2.25e12, 0.12)

function make_polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff; temp = 300.0, efield_freq = 0.0, verbose = true)

    # Collect data.
    ω = 2 * π * phonon_freq
    Ω = efield_freq
    T = temp
    α = frohlich_α(ϵ_optic, ϵ_static, phonon_freq, m_eff)

    # Prepare empty arrays for different temperatures.
    β = [] # Reduced thermodynamic beta (unitless)
    v = [] # Variational parameter v (1 / s, Hz)
    w = [] # Variational parameter w (1 / s, Hz)
    κ = [] # Spring constant (kg / s^2)
    M = [] # Fictitious mass (kg)
    F = [] # Free energy (meV)
    μ = [] # Mobility (cm^2 / Vs)
    Γ = [] # Optical absorption (1 / cm)

    # Initialise variation parameters.
    v_t, w_t = 0.0, 0.0
    if verbose
        print("\n")
    end

    for t in T # Iterate over temperatures.
        if verbose
            print("\e[2K", "Working on temperature: $(t) K / $(temp[end]) K.\n")
        end

        if t == 0.0 # If T = 0
            append!(β, Inf)  # set β = Inf

            # Evaluate variational parameters.
            v_t, w_t = variation(α; v = v_t, w = w_t)
            append!(v, v_t)
            append!(w, w_t)

            # Evaluate fictitious particle mass and spring constant.
            κ_t = (v_t^2 - w_t^2) # spring constant
            M_t = ((v_t^2 - w_t^2) / w_t^2) # mass
            append!(κ, κ_t)
            append!(M, M_t)

            # Evaluate free energy at zero temperature. NB: Enthalpy.
            F_t = free_energy(v_t, w_t, α) * 1000 * ħ * ω / eV
            append!(F, F_t)

            # Broadcast data.
            if verbose
                println("\e[2K", "α: ", round(α, digits = 3), " | β: ", Inf, " | v: ", round(v_t, digits = 3), " s^-1 | w: ", round(w_t, digits = 3), " s^-1")
                println("\e[2K", "κ: ", round(κ_t, digits = 3), " m_e kg/s^2 | M: ", round(M_t, digits = 3), " m_e kg")
                println("\e[2K", "Free Energy (Enthalpy): ", round(F_t, digits = 3), " meV")
            end

            # Prepare empty arrays for different frequencies.
            μ_t = [] # Mobilities
            Γ_t = [] # Optical absorptions

            for f in Ω # Iterate over frequencies

                if f < 1.0 # If Ω < 1 at T = 0

                    # Evaluate AC mobilities. NB: Ω = 0 is DC mobility.
                    μ_f = Inf # Infinite for Ω < 1 at T = 0
                    append!(μ_t, μ_f)

                    # Evaluate frequency-dependent optical absorptions.
                    Γ_f = 0.0 # Zero for Ω < 1 at T = 0
                    append!(Γ_t, Γ_f)

                    # Broadcast data.
                    if verbose
                        println("\e[2K", "Working on Frequency: $(f) Hz / $(efield_freq[end]) Hz")
                        println("\e[2K", "DC Mobility: ", round(μ_f, digits = 3), " cm^2/Vs")
                        println("\e[2K", "Optical absorption: ", round(Γ_f, digits = 3), " cm^-1")
                    end

                elseif f >= 1.0 # If Ω ≥ 1 at T = 0

                    # Evaluate AC mobilities. NB: β = 1000 for T ≈ 0.
                    μ_f = 100^2 * eV * polaron_mobility(f, 1e3, α, v_t, w_t) / (ω * m_eff * m_e)
                    append!(μ_t, μ_f)

                    # Evaluate optical absorptions. NB: β = 1000 for T ≈ 0. T = 0 unstable.
                    Γ_f = optical_absorption(f, 1e3, α, v_t, w_t) / (100 * c * ϵ_0 * sqrt(ϵ_optic))
                    append!(Γ_t, Γ_f)

                    # Broadcast data.
                    if verbose
                        println("\e[2K", "Working on Frequency: $(f) Hz / $(efield_freq[end]) Hz")
                        println("\e[2K", "AC Mobility: ", round(μ_f, digits = 3), " cm^2/Vs")
                        println("\e[2K", "Optical absorption: ", round(Γ_f, digits = 3), " cm^-1")
                    end
                end

                print("\033[F"^3) # Move cursor back
            end

            # Add freq-dependent mobilities and absorptions for T = 0.
            append!(μ, [μ_t])
            append!(Γ, [Γ_t])
            if verbose
                print("\033[F"^4) # Move cursor back
            end

        elseif t > 0.0 # If T > 0

            # Evaluate reduced thermodynamic beta.
            β_t = ħ * ω / (k_B * t)
            append!(β, β_t)

            # Evaluate variational parameters.
            v_t, w_t = variation(α, β_t; v = v_t, w = w_t)
            append!(v, v_t)
            append!(w, w_t)

            # Evaluate fictitious particle mass and spring constant.
            κ_t = (v_t^2 - w_t^2) # spring constant
            M_t = ((v_t^2 - w_t^2) / w_t^2) # mass
            append!(κ, κ_t)
            append!(M, M_t)

            # Evaluate free energy at finite temperature.
            F_t = free_energy(v_t, w_t, α, β_t) * 1000 * ħ * ω / eV
            append!(F, F_t)

            # Broadcast data.
            if verbose
                println("\e[2K", "α: ", round(α, digits = 3), " | β: ", round(β_t, digits = 3), " | v: ", round(v_t, digits = 3), " s^-1 | w: ", round(w_t, digits = 3), " s^-1")
                println("\e[2K", "κ: ", round(κ_t, digits = 3), " m_e kg/s^2 | M: ", round(M_t, digits = 3), " m_e kg")
                println("\e[2K", "Free Energy: ", round(F_t, digits = 3), " meV")
            end

            # Prepare empty arrays for different frequencies.
            μ_t = [] # Mobilities
            Γ_t = [] # Optical absorptions

            for f in Ω # Iterate over frequencies

                if f == 0.0 # If Ω = 0 at T > 0

                    # Evaluate DC mobility. NB: Ω = 0.001 ≈ 0. Ω = 0 unstable.
                    μ_f = 100^2 * eV * polaron_mobility(1e-3, β_t, α, v_t, w_t) / (ω * m_eff * m_e)
                    append!(μ_t, μ_f)

                    # Evaluate DC optical absorption. NB: Ω = 0.001 ≈ 0. Ω = 0 unstable.
                    Γ_f = optical_absorption(1e-3, β_t, α, v_t, w_t) / (100 * c * ϵ_0 * sqrt(ϵ_optic))
                    append!(Γ_t, Γ_f)

                    # Broadcast data.
                    if verbose
                        println("\e[2K", "Working on Frequency: $(f) Hz / $(efield_freq[end]) Hz")
                        println("\e[2K", "DC Mobility: ", round(μ_f, digits = 3), " cm^2/Vs")
                        println("\e[2K", "Optical absorption: ", round(Γ_f, digits = 3), " cm^-1")
                    end

                elseif f > 0.0 # If Ω > 0 at T > 0

                    # Evaluate AC mobilities.
                    μ_f = 100^2 * eV * polaron_mobility(f, β_t, α, v_t, w_t) / (ω * m_eff * m_e)
                    append!(μ_t, μ_f)

                    # Evaluate optical absorptions.
                    Γ_f = optical_absorption(f, β_t, α, v_t, w_t) / (100 * c * ϵ_0 * sqrt(ϵ_optic))
                    append!(Γ_t, Γ_f)

                    # Broadcast data.
                    if verbose
                        println("\e[2K", "Working on Frequency: $(f) Hz / $(efield_freq[end]) Hz")
                        println("\e[2K", "DC Mobility: ", round(μ_f, digits = 3), " cm^2/Vs")
                        println("\e[2K", "Optical absorption: ", round(Γ_f, digits = 3), " cm^-1")
                    end
                end

                if verbose
                    print("\033[F"^3) # Move cursor back
                end
            end

            # Add freq-dependent mobilities and absorptions for T > 0.
            append!(μ, [μ_t])
            append!(Γ, [Γ_t])
            if verbose
                print("\033[F"^4) # Move cursor back
            end

        else # Negative temperatures are unphysical!
            println("Temperature must be either zero or positive.")
        end
    end

    if verbose
        print("\n"^7) # Clear prints
    end

    # Return Polaron mutable struct with evaluated data.
    return Polaron(α, T, β, v, w, κ, M, F, Ω, μ, Γ)
end
