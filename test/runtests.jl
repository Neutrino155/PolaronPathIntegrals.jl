
using PolaronPathIntegrals
using Test

@testset "PolaronPathIntegrals" begin
    include("test_coupling.jl")
    include("test_free_energy.jl")
    include("test_variation.jl")
end
