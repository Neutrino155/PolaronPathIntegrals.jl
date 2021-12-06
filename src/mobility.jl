# mobility.jl

"""
----------------------------------------------------------------------
The Moblity of the Polaron.
----------------------------------------------------------------------
"""

"""
polaron_mobility(β::Float64, α::Float64, v::Float64, w::Float64; rtol = 1e-3)

    Calculate the dc mobility μ of the polaron at finite temperatues (equation (11.5) in [3]) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error for the integral to reach.
"""

function polaron_mobility(β, α, v, w; rtol = 1e-3)
	return 1 / imag(polaron_memory_function_dc(β, α, v, w; rtol = rtol))
end

"""
multi_mobility(β::Float64, α::Float64, v::Float64, w::Float64; rtol = 1e-3)

    Calculate the dc mobility μ of the polaron at finite temperatues (equation (11.5) in [3]) for a given frequency Ω. β is the thermodynamic beta. v and w are the variational polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. rtol specifies the relative error for the integral to reach.
"""

# function polaron_mobility(β, α, v, w; rtol = 1e-3)
# 	return 1 / imag(polaron_memory_function_dc(β, α, v, w; rtol = rtol))
# end