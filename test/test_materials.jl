

@testset "Materials" begin

S = Sphere()
R = GRAYTR.Ray(Point3(-1, 1/sqrt(2), 0), Vector3(1, 0, 0))
dg, a, b = GRAYTR.shape_intersect(R, S)
@test dg != nothing

@testset "World to local transformation" begin
    T = GRAYTR.local_transformation(dg)
    @test typeof(T) == Transformation
    
    invT = inv(T)
    
    @test T(dg.n) ≈ Normal3(0.0, 0.0, 1.0)
    @test invT(Normal3(0.0, 0.0, 1.0)) ≈ dg.n

    u = invT(Normal3(0.0, 1.0, 0.0))
    v = invT(Normal3(1.0, 0.0, 0.0))
    
    @test dot(u, dg.n) ≈ 0.0
    @test dot(v, dg.n) ≈ 0.0
    @test dot(u, v) ≈ 0.0
end


end # testset Materials
