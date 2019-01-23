
@testset "Dumb aggregate" begin


    sph1 = GRAYTR.Sphere()
    sph2 = GRAYTR.Sphere(translation(0, -2, 0))
    
    mat = GRAYTR.Lambert(GRAYTR.SingleLine(532.0, 1.0))
    
    P1 = GRAYTR.GeometricPrimitive(sph1, mat, nothing, 0)
    P2 = GRAYTR.GeometricPrimitive(sph2, mat, nothing, 0)
    
    Agg = GRAYTR.DumbAggregate([P1, P2])
    
    @test typeof(Agg) <: GRAYTR.DumbAggregate

    @test GRAYTR.world_bounds(Agg) == GRAYTR.BoundingBox(Point3(-1,-3,-1), Point3(1,1,1))
    
    R = Ray(Point3(5, 0, 0), Vector3(-1, 0, 0))
    @test GRAYTR.intersectP(R, Agg)
    isect = GRAYTR.intersect(R, Agg)
    @test isect.tmin ≈ 4.0
    @test isect.geometry.p ≈ Point3(1, 0, 0)
    @test isect.geometry.n ≈ Normal3(1, 0, 0)


    R = Ray(Point3(5, 2, 0), Vector3(-1, 0, 0))
    @test !GRAYTR.intersectP(R, Agg)
    isect = GRAYTR.intersect(R, Agg)
    @test isect == nothing


end
