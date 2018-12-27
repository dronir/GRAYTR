
import GRAYTR: lerp, quadratic

@testset "Utility functions" begin

@testset "Linear extrapolation/interpolation" begin
    @test lerp(0, 0, 1) == 0
    @test lerp(1, 0, 1) == 1
    @test lerp(0.0, 0, 1) ≈ 0.0
    @test lerp(1.0, 0, 1) ≈ 1.0

    for t = [-0.2, 0.15, 0.5, 1.25]
        @test lerp(t, 0.0, 1.0) ≈ t
        @test lerp(t, 0.0, 2.0) ≈ 2.0 * t
        @test lerp(t, -1.0, 1.0) ≈ 2.0 * t - 1.0
    end
end

@testset "Quadratic solver" begin
    @testset "No solution" begin
        solvable, t1, t2 = quadratic(2.0, 1.0, 1.0)
        @test !solvable
        @test t1 == t2 == Inf
    end
    
    @testset "Golden ratio" begin
        solvable, t1, t2 = quadratic(1.0, -1.0, -1.0)
        @test solvable
        @test t1 ≈ (1 - sqrt(5)) / 2
        @test t2 ≈ (1 + sqrt(5)) / 2
    end
end

end # outer testset
