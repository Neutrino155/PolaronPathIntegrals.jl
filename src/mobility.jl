# mobility.jl

"""
----------------------------------------------------------------------
The Moblity of the Polaron.
----------------------------------------------------------------------
"""

"""
polaron_mobility(Ω::Float64, β::Float64, α::Float64, v::Float64, w::Float64)

    Calculate the ac moblity μ(Ω) of the polaron at finite temperature (equation (1) in Hellwarth 1999) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling.
"""
function polaron_mobility(Ω, β, α, v, w)
    Ω / imag(χ(Ω, β, α, v, w))
end
