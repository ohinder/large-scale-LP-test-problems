include("utils.jl")
using Test
@testset "test rescale_instance" begin
    model = Model()
    @variable(model, x >= -2.0)
    @variable(model, -3.0 <= y <= 3)
    @objective(model, Min, 12x + 20y)
    @constraint(model, c1, 1.0x + 5.0y == 200.0)
    @constraint(model, c2, 6x + 8y >= 0.0)
    @constraint(model, c3, 7x + 12y >= 120)
    @constraint(model, c4, -2.0 <= 2x + 1y <= 0.0)

    lp_data = lp_matrix_data(model)
    rescaled_model = rescale_instance(lp_data)

    using SparseArrays, LinearAlgebra
    rescaled_data = lp_matrix_data(rescaled_model)

    @test rescaled_data.A ≈ sparse([ 
            [1.0 / (200 * 12)      5.0 / (20 * 200)]
            [6.0 / (12)    8 / (20.0)]
            [7.0 / (120.0 * 12)      12.0 / (20 * 120)]
            [2.0 / (12 * 2.0)      1.0 / (20 * 2.0)]
        ]) atol=1e-16
    @test rescaled_data.b_lower ≈ [1.0, 0.0, 1.0, -1.0] atol=1e-16
    @test rescaled_data.b_upper ≈ [1.0, Inf, Inf, 0.0]  atol=1e-16
    @test rescaled_data.c ≈ [1.0, 1.0] atol=1e-16
    @test rescaled_data.x_lower ≈ [-2.0 * 12, -3.0 * 20] atol=1e-16
    @test rescaled_data.x_upper ≈ [Inf, 3.0 * 20] atol=1e-16
end