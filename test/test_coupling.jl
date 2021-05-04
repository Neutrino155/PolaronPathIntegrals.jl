
@testset "frohlich_α" begin

    # NaCl ∼ 5 (Feynman 1955)
    @test frohlich_α(2.3, 5.6, 4.9e13 / (2 * π), 1.0) ≈ 5.0 atol = 0.3

    # CdTe ∼ 0.39 (Stone), 0.29 (Devreese)
    @test frohlich_α(7.1, 10.4, 5.08e12, 0.095) ≈ 0.3 atol = 0.1

    # GaAs ∼ 0.068 (Devreese)
    @test frohlich_α(10.89, 12.9, 8.46e12, 0.063) ≈ 0.068 atol = 0.01

end
