
@testset "Pressure camera" begin

@testset "Parallel ray generator" begin
    ray = GRAYTR.ray_parallel(2.0, 0.0, 0.0)
    @test ray.origin ≈ Point3(0.0, 0.0, 1.0)
    @test ray.direction ≈ Vector3(0.0, 0.0, -1.0)
    @test norm(ray.direction) ≈ 1.0
    
    # Need to explicitly loosen test tolerance here because of float instability
    # with trinometry. There are some vector components that should be zero but get
    # values on the order of 1e-16.
        
    ray = GRAYTR.ray_parallel(2.0, 1.0, 0.0)
    @test ray.direction ≈ Vector3(0.0, 0.0, -1.0)
    @test isapprox(ray.origin, Point3(2.0, 0.0, 1.0) ; atol=1e-15)
    
    ray = GRAYTR.ray_parallel(2.0, 0.0, 1.0)
    @test ray.direction ≈ Vector3(0.0, 0.0, -1.0)
    @test ray.origin ≈ Point3(0.0, 0.0, 1.0)
    
    ray = GRAYTR.ray_parallel(2.0, 1.0, 1.0)
    @test ray.direction ≈ Vector3(0.0, 0.0, -1.0)
    @test isapprox(ray.origin, Point3(2.0, 0.0, 1.0) ; atol=1e-15)
    
    ray = GRAYTR.ray_parallel(2.0, 1.0, 0.5)
    @test isapprox(ray.origin, Point3(-2.0, 0.0, 1.0) ; atol=1e-15)

end

@testset "Point-source ray generator" begin
    ray = GRAYTR.ray_from_point(2.0, 4.0, 1.0, 0.0)
    @test ray.origin ≈ Point3(0.0, 0.0, 4.0)
    @test ray.direction ≈ Vector3(0.0, 0.0, -1.0)
    @test norm(ray.direction) ≈ 1.0
    
    # Need to explicitly loosen test tolerance here because of float instability
    # with trinometry. There are some vector components that should be zero but get
    # values on the order of 1e-16.
    
    ray = GRAYTR.ray_from_point(2.0, 4.0, 1.0, 1.0)
    @test ray.origin ≈ Point3(0.0, 0.0, 4.0)
    @test ray.direction ≈ Vector3(0.0, 0.0, -1.0)
    @test norm(ray.direction) ≈ 1.0

    ct = sqrt(1.0 - (2.0 / 4.0)^2)
    
    ray = GRAYTR.ray_from_point(2.0, 4.0, 0.0, 0.0)
    @test ray.origin ≈ Point3(0.0, 0.0, 4.0)
    @test isapprox(ray.direction, Vector3(sqrt(1-ct^2), 0.0, -ct), atol=1e-15)
    @test norm(ray.direction) ≈ 1.0
    
    ray = GRAYTR.ray_from_point(2.0, 4.0, 0.0, 0.5)
    @test ray.origin ≈ Point3(0.0, 0.0, 4.0)
    @test isapprox(ray.direction, Vector3(-sqrt(1-ct^2), 0.0, -ct), atol=1e-15)
    @test norm(ray.direction) ≈ 1.0
    
    ray = GRAYTR.ray_from_point(2.0, 4.0, 0.0, 0.25)
    @test ray.origin ≈ Point3(0.0, 0.0, 4.0)
    @test isapprox(ray.direction, Vector3(0.0, sqrt(1-ct^2), -ct), atol=1e-15)
    @test norm(ray.direction) ≈ 1.0
    
    ray = GRAYTR.ray_from_point(2.0, 4.0, 0.0, 0.75)
    @test ray.origin ≈ Point3(0.0, 0.0, 4.0)
    @test isapprox(ray.direction, Vector3(0.0, -sqrt(1-ct^2), -ct), atol=1e-15)
    @test norm(ray.direction) ≈ 1.0
    

end

end

