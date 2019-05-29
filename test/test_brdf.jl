@testset "BRDF" begin

@testset "Utilities" begin
    w = normalize(Vector3(1, 0, 1))
    
    @test GRAYTR.costheta(w) ≈ sqrt(2)/2
    @test GRAYTR.sintheta(w) ≈ sqrt(2)/2
    @test GRAYTR.sintheta2(w) ≈ 0.5
    @test GRAYTR.cosphi(w) ≈ 1.0
    @test GRAYTR.sinphi(w) ≈ 0.0
end



@testset "Materials" begin

S = Sphere()
R = GRAYTR.Ray(Point3(-1, 1/sqrt(2), 0), Vector3(1, 0, 0))
dg, a, b = GRAYTR.shape_intersect(R, S)
@test dg != nothing

@testset "World to local transformation" begin
    @test typeof(dg.T) == Transformation
    
    invT = inv(dg.T)
    
    @test dg.T(dg.n) ≈ Normal3(0.0, 0.0, 1.0)
    @test invT(Normal3(0.0, 0.0, 1.0)) ≈ dg.n

    u = invT(Normal3(0.0, 1.0, 0.0))
    v = invT(Normal3(1.0, 0.0, 0.0))
    
    @test dot(u, dg.n) ≈ 0.0
    @test dot(v, dg.n) ≈ 0.0
    @test dot(u, v) ≈ 0.0
end


end # testset Materials


spec = GRAYTR.SingleLine(532.0, 1.0)
w0 = normalize(Vector3(0, 0, 1))
w1 = normalize(Vector3(1, 0, 1))
ws = normalize(Vector3(-1, 0, 1))

@testset "Lambert" begin
    R = GRAYTR.Lambert(spec)
    
    @test GRAYTR.evaluate(R, w0, w1) == GRAYTR.SingleLine(532.0, 1) / π
    @test GRAYTR.rho(R, w0) == GRAYTR.SingleLine(532.0, 1)
end

@testset "Lommel-Seeliger" begin
    P(alpha) = 1.0
    R = GRAYTR.LommelSeeliger(spec, P)
    
    @test GRAYTR.evaluate(R, w0, w1) == GRAYTR.SingleLine(532.0, 1) * 0.25 / (1.0 + sqrt(2)/2)
    
    
end

@testset "Ashkhmin-Shirley" begin
#    R = GRAYTR.AshkhminShirleySingle(spec, 0.5, 100.0)
#    @test GRAYTR.evaluate(R, w0, w1) <: GRAYTR.Spectrum    
end


@testset "SpecularDiffuse" begin 
    L = GRAYTR.Lambert(spec)
    R = GRAYTR.SpecularDiffuse(spec, 1.0, 0.0)
    
    @test GRAYTR.evaluate(R, w0, w1) == GRAYTR.SingleLine(532.0, 1) / π
    @test GRAYTR.evaluate(R, w0, w1) == GRAYTR.evaluate(L, w0, w1)
    
    R = GRAYTR.SpecularDiffuse(spec, 0.5, 0.5)
    @test GRAYTR.evaluate(R, w1, ws) == GRAYTR.SingleLine(532.0, 1) * (0.5 + 0.5/π)
    
    @test GRAYTR.compute_pressure(R, w0, spec) ≈ Vector3(0.0, 0.0, -(1 + 5/6))
end



end # main BRDF testset

