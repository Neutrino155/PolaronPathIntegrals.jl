# free_energy.jl

"""
Implementation of Feynman's original variational technique applied to the Polaron model.

See Feynman 1955:
http://dx.doi.org/10.1103/PhysRev.97.660
"""

# Equation 31: The <|X(t) - X(s)|^{-1}> * exp(-|t-w|) effective action.
A_integrand(x, v, w) = (w^2 * x + (v^2 - w^2) / v * (1 - exp(-v * x)))^(-0.5) * exp(-x)

A(v, w, α) = π^(-0.5) * α * v * QuadGK.quadgk(x -> A_integrand(x, v, w), 0, Inf)[1]

# Equation 33: Lowest Free energy E = -B - A where B = -3/(4v)*(v-w)^2.
free_energy(v, w, α) = (3 / (4 * v)) * (v - w)^2 - A(v, w, α)

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
A(v, w, α, β) = α * v / (sqrt(π) * (exp(BigFloat(β)) - 1)) * QuadGK.quadgk(x -> A_integrand(x, v, w, β), BigFloat(0), BigFloat(β / 2))[1]

# Equation 62b in Hellwarth. Equation 20 in Osaka.
B(v, w, β) = 3 / β * (log(v / w) - 1 / 2 * log(2 * π * BigFloat(β)) - log(sinh(v * BigFloat(β) / 2) / sinh(w * BigFloat(β) / 2)))

# Equation 62e in Hellwarth. Equation 17 in Osaka.
C(v, w, β) = 3 / 4 * (v^2 - w^2) / v * (coth(v * BigFloat(β) / 2) - 2 / (v * β))

# Equation 62a in Hellwarth. In paragraph below Equation 22 in Osaka; has extra 1/β due to different definition of A, B & C.
function free_energy(v, w, α, β)
    setprecision(BigFloat, 64)
    a = A(v, w, α, β)
    b = B(v, w, β)
    c = C(v, w, β)
    -(a + b + c)
end

"""
Multiple Phonon Branches
"""

function κ_i(i, v, w) # fictitious spring constant, multiple variational params
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

function h_i(i, v, w) # some vector relating to harmonic eigenmodes
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

function C_ij(i, j, v, w) # generalised Feynman C variational parameter (inclusive of multiple v and w params)
    C = w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
    return C
end

function D_j(τ, β, v, w) # log of dynamic structure factor for polaron 
    D = τ * (1 - τ / β)
    for i in 1:length(v)
        if v[i] != w[i]
        D += (h_i(i, v, w) / v[i]^2) * (2 * sinh(v[i] * τ / 2) * sinh(v[i] * (β - τ) / 2) / (v[i] * sinh(v[i] * β / 2)) - τ * (1 - τ / β))
        end
    end
    return D
end

function multi_free_energy(v_params, w_params, T, ϵ_optic, m_eff, volume, freqs_and_ir_activity, phonon_branch)

    setprecision(BigFloat, 32) # Speed up. Stops potential overflows.

    # Extract phonon frequencies and ir activities.
    phonon_freqs = freqs_and_ir_activity[:, 1]
    ir_activity = freqs_and_ir_activity[:, 2]

    num_of_branches = length(phonon_freqs)
    j = phonon_branch # jth phonon branch
    
    # total dielectric contribution from all phonon branches (used as a normalisation)
    ϵ_tot = ϵ_total(freqs_and_ir_activity, volume)

    # Generalisation of B i.e. Equation 62c in Hellwarth.
    S_integrand(τ, β, v, w) = cosh(β / 2 - τ) / (sinh(β / 2) * sqrt(abs(D_j(τ, β, v, w))))
    S(β, α, v, w) = α / √π * quadgk(τ -> S_integrand(τ, β, v, w), 0.0, β / 2)[1]

    # Generalisation of C i.e. Equation 62e in Hellwarth.
    function S_0(β, v, w)
        s = 0.0
        for i in 1:length(v)
            for j in 1:length(v)
                s += C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2)  - 2 / (β * v[j]))
            end
        end
        3 * s / num_of_branches 
    end

    # Generalisation of A i.e. Equation 62b in Hellwarth.
    function log_Z_0(β, v, w)
        s = -log(2π * β) / 2
        for i in 1:length(v)
            if v[i] != w[i]
                s += log(v[i] / w[i]) -log(sinh(v[i] * β / 2) / sinh(w[i] * β / 2))
            end
        end
        3 / β * s / num_of_branches
    end

    ω = 2π * 1e12 * phonon_freqs[j] # angular phonon freq im 2π Hz
    β = BigFloat(ħ * ω / (kB * T)) # reduced thermodynamic beta
    ϵ_ionic = ϵ_ionic_mode(phonon_freqs[j], ir_activity[j], volume) # ionic dielectric contribution for current phonon branch
    α = frohlich_α_j(ϵ_optic, ϵ_ionic, ϵ_tot, phonon_freqs[j], m_eff) # decomposed alpha for current phonon branch

    # Print out data.
    # println("α = $α, β = $β, f = $(phonon_freqs[j]), ")
    # print("F0 = $(log_Z_0(β, ω, v_params .* ω, w_params .* ω) / β), ")
    # print("<S0> = $(S_0(β, ω, v_params .* ω, w_params .* ω) / β), ")
    # print("<S> = $(S(β, α, ω, v_params .* ω, w_params .* ω) / β), ")
    # println("F = $(-log_Z_0(β, ω, v_params .* ω, w_params .* ω) / β - S(β, α, ω, v_params .* ω, w_params .* ω) / β + S_0(β, ω, v_params .* ω, w_params .* ω) / β)")
    
    # F = -(A + B + C) in Hellwarth.
    F = -(S(β, α, v_params, w_params) + S_0(β, v_params, w_params) + log_Z_0(β, v_params, w_params)) * ħ * ω # × ħω branch phonon energy
    
    return F / eV * 1e3 # change to meV
end