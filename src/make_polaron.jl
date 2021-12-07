# make_polaron.jl

# Example: MAPI: make_polaron(4.5, 24.1, 2.25e12, 0.12)

function make_polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff; temp = 300.0, efield_freq = 0.0, volume = nothing, ir_activity = nothing, N_params = 1, verbose = true)

    # Collect data.
    Ω = efield_freq
    T = temp
    N_modes = length(phonon_freq)
    N_temp = length(T)
    N_freq = length(efield_freq)
    
    # If ir_activity is given (an array) then multiple phonon modes are present.
    if N_modes == 1
        ω = 2π * phonon_freq
        α = frohlich_α(ϵ_optic, ϵ_static, phonon_freq, m_eff)
    else
        ω = 2π .* phonon_freq
        ϵ_ionic = [ϵ_ionic_mode(i, j, volume) for (i, j) in zip(phonon_freq, ir_activity)]
        α = [frohlich_α(ϵ_optic, i, sum(ϵ_ionic), j, m_eff) for (i, j) in zip(ϵ_ionic, phonon_freq)]
    end

    # Prepare empty arrays for different temperatures.
    β = Matrix{Float64}(undef, N_temp, N_modes) # Reduced thermodynamic beta (unitless)
    v = Vector{Float64}(undef, N_temp) # Variational parameter v (1 / s, Hz)
    w = Vector{Float64}(undef, N_temp) # Variational parameter w (1 / s, Hz)
    κ = Vector{Float64}(undef, N_temp) # Spring constant (kg / s^2)
    M = Vector{Float64}(undef, N_temp) # Fictitious mass (kg)
    F = Vector{Float64}(undef, N_temp) # Free energy (meV)
    Z = Matrix{ComplexF64}(undef, N_freq, N_temp) # Complex impedence (cm^2 / Vs)
    σ = Matrix{ComplexF64}(undef, N_freq, N_temp) # Complex conductivity (1 / cm)

    # Initialise variation parameters.
    v_t, w_t = 0.0, 0.0
    if verbose
        print("\n")
    end

    for t in 1:length(T) # Iterate over temperatures.
        if verbose
            print("\e[2K", "Working on temperature: $(T[t]) K / $(T[end]) K.\n")
        end

        if T[t] == 0.0 # If T = 0
            β[t, :] = repeat([Inf], N_modes)  # set β = Inf

            # Evaluate variational parameters.
            v_t, w_t = variation(α; v = v_t, w = w_t, ω = ω)

            v[t] = v_t
            w[t] = w_t

            # Evaluate fictitious particle mass and spring constant.
            κ_t = (v_t^2 - w_t^2) # spring constant
            M_t = ((v_t^2 - w_t^2) / w_t^2) # mass
            κ[t] = κ_t
            M[t] = M_t

            # Evaluate free energy at zero temperature. NB: Enthalpy.
            F_t = Float64(free_energy(v_t, w_t, α; ω = ω)) * 1000 * ħ / eV
            F[t] = F_t

            # Broadcast data.
            if verbose
                println("\e[2K", "α: ", round(sum(α), digits = 3), " | β: ", Inf, " | v: ", round(v_t, digits = 3), " s^-1 | w: ", round(w_t, digits = 3), " s^-1")
                println("\e[2K", "κ: ", round(κ_t, digits = 3), " m_e kg/s^2 | M: ", round(M_t, digits = 3), " m_e kg")
                println("\e[2K", "Free Energy (Enthalpy): ", round(F_t, digits = 3), " meV")
            end

            for f in 1:length(Ω) # Iterate over frequencies

                if Ω[f] == 0.0 # If Ω = 0 at T = 0

                    # Evaluate AC mobilities. NB: Ω = 0 is DC mobility.
                    Z_f = 0 + 1im * 0
                    Z[f, t] = Z_f

                    # Evaluate frequency-dependent optical absorptions.
                    σ_f = 0 + 1im * 0
                    σ[f, t] = σ_f

                    # Broadcast data.
                    if verbose
                        println("\e[2K", "Working on Frequency: $(Ω[f]) Hz / $(Ω[end]) Hz")
                        println("\e[2K", "DC Impedence: ", round(Z_f, digits = 3), " cm^2/Vs")
                        println("\e[2K", "DC Conductivity: ", round(σ_f, digits = 3), " cm^-1")
                    end

                elseif Ω[f] > 0.0 # If Ω > 0 at T = 0

                    # Evaluate AC mobilities.
                    Z_f = polaron_complex_impedence(Ω[f], β[t], α, v_t, w_t; ω = ω) / eV^2 * (m_e * m_eff) * 100^2 / 1e12
                    Z[f, t] = Z_f

                    # Evaluate optical absorptions.
                    σ_f = polaron_complex_conductivity(Ω[f], β[t], α, v_t, w_t; ω = ω) * eV^2 / (m_e * m_eff) / 100^2 * 1e12
                    σ[f, t] = σ_f

                    # Broadcast data.
                    if verbose
                        println("\e[2K", "Working on Frequency: $(Ω[f]) Hz / $(Ω[end]) Hz")
                        println("\e[2K", "AC Impedence: ", round(Z_f, digits = 3), " cm^2/Vs")
                        println("\e[2K", "AC Conductivity: ", round(σ_f, digits = 3), " cm^-1")
                    end
                end

                print("\033[F"^3) # Move cursor back
            end

            if verbose
                print("\033[F"^4) # Move cursor back
            end

        elseif T[t] > 0.0 # If T > 0

            # Evaluate reduced thermodynamic beta.
            β_t = ω .* ħ / (k_B * T[t]) * 1e12
            β[t, :] .= β_t

            # Evaluate variational parameters.
            v_t, w_t = variation(α, β_t; v = v_t, w = w_t, ω = ω)
            v[t] = v_t
            w[t] = w_t

            # Evaluate fictitious particle mass and spring constant.
            κ_t = (v_t^2 - w_t^2) # spring constant
            M_t = ((v_t^2 - w_t^2) / w_t^2) # mass
            κ[t] = κ_t
            M[t] = M_t

            # Evaluate free energy at finite temperature.
            F_t = Float64(free_energy(v_t, w_t, α, β_t; ω = ω)) * 1000 * ħ / eV
            F[t] = F_t

            # Broadcast data.
            if verbose
                println("\e[2K", "α: ", round(sum(α), digits = 3), " | β: ", round.(β_t, digits = 3), " | v: ", round(v_t, digits = 3), " s^-1 | w: ", round(w_t, digits = 3), " s^-1")
                println("\e[2K", "κ: ", round(κ_t, digits = 3), " m_e kg/s^2 | M: ", round(M_t, digits = 3), " m_e kg")
                println("\e[2K", "Free Energy: ", round(F_t, digits = 3), " meV")
            end

            for f in 1:length(Ω) # Iterate over frequencies

                if Ω[f] == 0.0 # If Ω = 0 at T > 0

                    # Evaluate DC mobility.
                    Z_f = polaron_complex_impedence(Ω[f], β[t, :], α, v_t, w_t; ω = ω) / eV^2 * (m_e * m_eff) * 100^2 / 1e12
                    Z[f, t] = Z_f

                    # Evaluate DC optical absorption. 
                    σ_f = polaron_complex_conductivity(Ω[f], β[t, :], α, v_t, w_t; ω = ω) * eV^2 / (m_e * m_eff) / 100^2 * 1e12
                    σ[f, t] = σ_f

                    # Broadcast data.
                    if verbose
                        println("\e[2K", "Working on Frequency: $(Ω[f]) Hz / $(Ω[end]) Hz")
                        println("\e[2K", "DC Impedence: ", round(Z_f, digits = 3), " cm^2/Vs")
                        println("\e[2K", "DC Conductivity: ", round(σ_f, digits = 3), " cm^-1")
                    end

                elseif Ω[f] > 0.0 # If Ω > 0 at T > 0

                    # Evaluate AC mobilities.
                    Z_f = polaron_complex_impedence(Ω[f], β[t, :], α, v_t, w_t; ω = ω) / eV^2 * (m_e * m_eff) * 100^2 / 1e12
                    Z[f, t] = Z_f

                    # Evaluate optical absorptions.
                    σ_f = polaron_complex_conductivity(Ω[f], β[t, :], α, v_t, w_t; ω = ω) * eV^2 / (m_e * m_eff) / 100^2 * 1e12
                    σ[f, t] = σ_f

                    # Broadcast data.
                    if verbose
                        println("\e[2K", "Working on Frequency: $(Ω[f]) Hz / $(Ω[end]) Hz")
                        println("\e[2K", "AC Impedence: ", round(Z_f, digits = 3), " cm^2/Vs")
                        println("\e[2K", "AC Conductivity: ", round(σ_f, digits = 3), " cm^-1")
                    end
                end

                if verbose
                    print("\033[F"^3) # Move cursor back
                end
            end

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
    return Polaron(sum(α), T, β, v, w, κ, M, F, Ω, Z, σ)
end

function make_polaron(α; temp = 300.0, efield_freq = 0.0, verbose = true)

# Collect data.
Ω = efield_freq
T = temp
ω = 1

# Prepare empty arrays for different temperatures.
β = Vector{Float64}(undef, length(T)) # Reduced thermodynamic beta (unitless)
v = Vector{Float64}(undef, length(T)) # Variational parameter v (1 / s, Hz)
w = Vector{Float64}(undef, length(T)) # Variational parameter w (1 / s, Hz)
κ = Vector{Float64}(undef, length(T)) # Spring constant (kg / s^2)
M = Vector{Float64}(undef, length(T)) # Fictitious mass (kg)
F = Vector{Float64}(undef, length(T)) # Free energy (meV)
Z = Matrix{ComplexF64}(undef, length(Ω), length(T)) # Mobility (cm^2 / Vs)
σ = Matrix{ComplexF64}(undef, length(Ω), length(T)) # Optical absorption (1 / cm)

# Initialise variation parameters.
v_t, w_t = 0.0, 0.0
if verbose
    print("\n")
end

for t in 1:length(T) # Iterate over temperatures.
    if verbose
        print("\e[2K", "Working on temperature: $(T[t]) K / $(T[end]) K.\n")
    end

    if T[t] == 0.0 # If T = 0
        β[t] = Inf  # set β = Inf

        # Evaluate variational parameters.
        v_t, w_t = variation(α; v = v_t, w = w_t)
        v[t] = v_t
        w[t] = w_t

        # Evaluate fictitious particle mass and spring constant.
        κ_t = (v_t^2 - w_t^2) # spring constant
        M_t = ((v_t^2 - w_t^2) / w_t^2) # mass
        κ[t] = κ_t
        M[t] = M_t

        # Evaluate free energy at zero temperature. NB: Enthalpy.
        F_t = Float64(free_energy(v_t, w_t, α))
        F[t] = F_t

        # Broadcast data.
        if verbose
            println("\e[2K", "α: ", round(α, digits = 3), " | β: ", Inf, " | v: ", round(v_t, digits = 3), " | w: ", round(w_t, digits = 3))
            println("\e[2K", "κ: ", round(κ_t, digits = 3), " | M: ", round(M_t, digits = 3))
            println("\e[2K", "Free Energy (Enthalpy): ", round(F_t, digits = 3))
        end

        for f in 1:length(Ω)# Iterate over frequencies

            if Ω[f] == 0.0 # If Ω = 0 at T = 0

                # Evaluate DC mobilities. NB: Ω = 0 is DC mobility.
                Z_f = 0 + 1im * 0
                Z[f, t] = Z_f

                # Evaluate frequency-dependent optical absorptions.
                σ_f = 0 + 1im * 0
                σ[f, t] = σ_f

                # Broadcast data.
                if verbose
                    println("\e[2K", "Working on Frequency: $(Ω[f]) Hz / $(Ω[end]) Hz")
                    println("\e[2K", "DC Mobility: ", round(Z_f, digits = 3))
                    println("\e[2K", "Optical absorption: ", round(σ_f, digits = 3))
                end

            elseif Ω[f] > 0.0 # If Ω > 0 at T = 0

                # Evaluate AC mobilities.
                Z_f = polaron_complex_impedence(Ω[f], Inf, α, v_t, w_t) 
                Z[f, t] = Z_f

                # Evaluate optical absorptions.
                σ_f = complex_conductivity(Ω[f], α, v_t, w_t)
                σ[f, t] = σ_f

                # Broadcast data.
                if verbose
                    println("\e[2K", "Working on Frequency: $(Ω[f]) Hz / $(Ω[end]) Hz")
                    println("\e[2K", "AC Mobility: ", round(Z_f, digits = 3))
                    println("\e[2K", "Optical absorption: ", round(σ_f, digits = 3))
                end
            end

            print("\033[F"^3) # Move cursor back
        end

        if verbose
            print("\033[F"^4) # Move cursor back
        end

    elseif T[t] > 0.0 # If T > 0

        # Evaluate reduced thermodynamic beta.
        β_t = ω / T[t]
        β[t] = β_t

        # Evaluate variational parameters.
        v_t, w_t = variation(α, β_t; v = v_t, w = w_t)
        v[t] = v_t
        w[t] = w_t

        # Evaluate fictitious particle mass and spring constant.
        κ_t = (v_t^2 - w_t^2) # spring constant
        M_t = ((v_t^2 - w_t^2) / w_t^2) # mass
        κ[t] = κ_t
        M[t] = M_t

        # Evaluate free energy at finite temperature.
        F_t = Float64(free_energy(v_t, w_t, α, β_t))
        F[t] = F_t

        # Broadcast data.
        if verbose
            println("\e[2K", "α: ", round(α, digits = 3), " | β: ", round(β_t, digits = 3), " | v: ", round(v_t, digits = 3), " | w: ", round(w_t, digits = 3))
            println("\e[2K", "κ: ", round(κ_t, digits = 3), " | M: ", round(M_t, digits = 3))
            println("\e[2K", "Free Energy: ", round(F_t, digits = 3))
        end

        for f in 1:length(Ω) # Iterate over frequencies

            if Ω[f] == 0.0 # If Ω = 0 at T > 0

                # Evaluate DC mobility.
                Z_f = complex_impedence_dc(β_t, α, v_t, w_t)
                Z[f, t] = Z_f

                # Evaluate DC optical absorption. 
                σ_f = complex_conductivity_dc(β_t, α, v_t, w_t)
                σ[f, t] = σ_f

                # Broadcast data.
                if verbose
                    println("\e[2K", "Working on Frequency: $(Ω[f]) Hz / $(Ω[end]) Hz")
                    println("\e[2K", "DC Mobility: ", round(Z_f, digits = 3))
                    println("\e[2K", "Optical absorption: ", round(σ_f, digits = 3))
                end

            elseif Ω[f] > 0.0 # If Ω > 0 at T > 0

                # Evaluate AC mobilities.
                Z_f = complex_impedence(Ω[f], β_t, α, v_t, w_t)
                Z[f, t] = Z_f

                # Evaluate optical absorptions.
                σ_f = complex_conductivity(Ω[f], β_t, α, v_t, w_t)
                σ[f, t] = σ_f

                # Broadcast data.
                if verbose
                    println("\e[2K", "Working on Frequency: $(Ω[f]) Hz / $(Ω[end]) Hz")
                    println("\e[2K", "AC Mobility: ", round(Z_f, digits = 3))
                    println("\e[2K", "Optical absorption: ", round(σ_f, digits = 3))
                end
            end

            if verbose
                print("\033[F"^3) # Move cursor back
            end
        end

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
return Polaron(α, T, β, v, w, κ, M, F, Ω, Z, σ)
end

function save_polaron(polaron; path = "../data/")
    polaron_df = DataFrames.DataFrame(
        alpha = [polaron.α for i in 1:length(polaron.T)],
        temperature = polaron.T,
        beta = polaron.β,
        v = polaron.v,
        w = polaron.w,
        spring_constant = polaron.κ,
        polaron_mass = polaron.M,
        free_energy = polaron.F,
    )
    mobility_df = DataFrames.DataFrame([[0.0, polaron.T...]'; [polaron.Ω polaron.μ]], :auto)
    absorption_df = DataFrames.DataFrame([[0.0, polaron.T...]'; [polaron.Ω polaron.Γ]], :auto)
    CSV.write(path * "polaron_data.csv", polaron_df)
    CSV.write(path * "mobility_data.csv", mobility_df)
    CSV.write(path * "absorption_data.csv", absorption_df)

end

function load_polaron(polaron_data_path, mobility_data_path, absorption_data_path)
    polaron_data = CSV.File(polaron_data_path) |> Tables.matrix
    mobility_data = CSV.File(mobility_data_path) |> Tables.matrix
    absorption_data = CSV.File(absorption_data_path) |> Tables.matrix
    return Polaron(
        polaron_data[1], 
        polaron_data[:, 2],
        polaron_data[:, 3],
        polaron_data[:, 4],
        polaron_data[:, 5],
        polaron_data[:, 6],
        polaron_data[:, 7],
        polaron_data[:, 8],
        mobility_data[2:end, 1],
        mobility_data[2:end, 2:end],
        absorption_data[2:end, 2:end]
    )
end