# mobility.jl

"""
----------------------------------------------------------------------
The Moblity of the Polaron.
----------------------------------------------------------------------
"""

"""
polaron_mobility_ac(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the ac moblity μ(Ω) of the polaron at finite temperature (equation (1) in Hellwarth 1999) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function polaron_mobility_ac(Ω, β, α, v, w)
    Ω / imag(χ(Ω, β, α, v, w))
end

"""
polaron_mobility_ac(Ω::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the ac moblity μ(Ω) of the polaron at zero-temperatures (equation (46) in FHIP 1962) for a given frequency Ω. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function polaron_mobility_ac(Ω, α, v, w)
    Ω / imag(χ(Ω, α, v, w))
end

"""
polaron_mobility_dc(β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the dc moblity μ(β) of the polaron at finite temperatures (equation (46) in FHIP 1962) at zero frequenc. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function polaron_mobility_dc(β, α, v, w)
    1 / imag(χ_dc(β, α, v, w))
end
