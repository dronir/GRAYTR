
@testset "Filters" begin

@testset "Box filter" begin
    # Box filter always returns 1.0
    B = GRAYTR.BoxFilter(1.0, 2.0)
    @test GRAYTR.evaluate(B, 0.0, 1.0) ≈ 1.0
    @test GRAYTR.evaluate(B, 2.0, -1.0) ≈ 1.0
end


@testset "Triangle filter" begin
    T = GRAYTR.TriangleFilter(1.0, 2.0)
    @test GRAYTR.evaluate(T, 0.0, 0.0) ≈ 2.0
    @test GRAYTR.evaluate(T, 0.5, 0.0) ≈ 1.0
    @test GRAYTR.evaluate(T, 0.5, 1.0) ≈ 0.5
    @test GRAYTR.evaluate(T, 10.0, 0.0) ≈ 0.0
end

end
