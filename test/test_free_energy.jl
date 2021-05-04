
athermal_values = [
    [3.00, 5.00, 7.00, 9.00, 11.00],
    [3.44, 4.02, 5.81, 9.85, 15.5],
    [2.55, 2.13, 1.6, 1.28, 1.15],
    [-3.1333, -5.4401, -8.1127, -11.486, -15.710]
]

@testset "free energy" begin

    @testset "athermal" begin
        for (α, v, w, E) in zip(athermal_values...)
            @test free_energy(v, w, α) ≈ E atol = 0.001
        end
    end

end
