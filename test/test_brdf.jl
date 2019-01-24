@testset "BRDF" begin

@testset "Utilities" begin
    w = normalize(Vector3(1, 0, 1))
    
    @test GRAYTR.costheta(w) ≈ sqrt(2)/2
    @test GRAYTR.sintheta(w) ≈ sqrt(2)/2
    @test GRAYTR.sintheta2(w) ≈ 0.5
    @test GRAYTR.cosphi(w) ≈ 1.0
    @test GRAYTR.sinphi(w) ≈ 0.0
end

@testset "Lambert" begin
    R = GRAYTR.Lambert(GRAYTR.SingleLine(532.0, 1.0))
    
    w0 = normalize(Vector3(0, 0, 1))
    w1 = normalize(Vector3(1, 0, 1))
    
    @test GRAYTR.BSDF_type(R) == 9
    @test GRAYTR.evaluate(R, w0, w1) == GRAYTR.SingleLine(532.0, 1) / π
    
    @test GRAYTR.rho(R, w0) == GRAYTR.SingleLine(532.0, 1)
end

@testset "Lommel-Seeliger" begin
end

@testset "Ashkhmin-Shirley" begin
end

end

