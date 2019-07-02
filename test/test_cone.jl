@testset "Cone" begin

cone = Cone()
@test cone.id == 1
@test cone.h == 1.0
@test GRAYTR.can_intersect(cone)

cone = Cone(2)
@test cone.id == 2
@test cone.h == 1.0

# Create with transformation
cone = Cone(GRAYTR.scaling(2.0, 1.0, 1.0))
@test GRAYTR.obj_bounds(cone) == GRAYTR.BoundingBox(Point3(-1, -1, 0), Point3(1, 1, 1))
@test GRAYTR.world_bounds(cone) == GRAYTR.BoundingBox(Point3(-2, -1, 0), Point3(2, 1, 1))

cone = Cone(2, GRAYTR.scaling(2.0, 1.0, 1.0))
@test cone.id == 2

# Transform after creation
T = GRAYTR.scaling(1.0, 2.0, 1.0)
cone = T(Cone())
@test GRAYTR.obj_bounds(cone) == GRAYTR.BoundingBox(Point3(-1, -1, 0), Point3(1, 1, 1))
@test GRAYTR.world_bounds(cone) == GRAYTR.BoundingBox(Point3(-1, -2, 0), Point3(1, 2, 1))



@testset "Intersections" begin

    # intersectP
    D = Cone(translation(0, 0, -1))
    
    R = GRAYTR.Ray(Point3(0.5, 0, 1), Vector3(0, 0, -1))
    @test GRAYTR.intersectP(R, D)
    
    R = GRAYTR.Ray(Point3(0.5, 0, 1), Vector3(0, 0, -1), 10, Inf, 1)
    @test !GRAYTR.intersectP(R, D)
    
    R = GRAYTR.Ray(Point3(0, 2, 1), Vector3(0, 0, -1))
    @test !GRAYTR.intersectP(R, D)
    
    R = GRAYTR.Ray(Point3(0.9*sqrt(2), 0.9*sqrt(2), 1), Vector3(0, 0, -1))
    @test !GRAYTR.intersectP(R, D)
    
    R = GRAYTR.Ray(Point3(0.5, 0, -2), Vector3(0, 0, 1))
    @test GRAYTR.intersectP(R, D)
    
    R = GRAYTR.Ray(Point3(0.0, 0, -0.2), Vector3(0, 1, 0))
    @test GRAYTR.intersectP(R, D)
    
    R = GRAYTR.Ray(Point3(2.0, 0, 0.2), Vector3(-1, 0, 0))
    @test !GRAYTR.intersectP(R, D)
    
    
    # intersect
    
    R = GRAYTR.Ray(Point3(0.5, 0, 1), Vector3(0, 0, -1))
    dg, t, e = GRAYTR.shape_intersect(R, D)
    @test dg.p ≈ Point3(0.5, 0, -0.5)
    @test dg.n ≈ normalize(Normal3(1, 0, 1))
    @test t ≈ 1.5
    
    R = GRAYTR.Ray(Point3(0.5, 0, -3), Vector3(0, 0, 1))
    dg, t, e = GRAYTR.shape_intersect(R, D)
    @test dg.p ≈ Point3(0.5, 0, -0.5)
    @test dg.n ≈ normalize(Normal3(-1.0, 0, -1.0))
    @test t ≈ 2.5
    
    R = GRAYTR.Ray(Point3(0.9*sqrt(2), 0.9*sqrt(2), 2), Vector3(0, 0, -1))
    dg, t, e = GRAYTR.shape_intersect(R, D)
    @test dg == nothing
    @test isnan(t)
    @test isnan(e)
    
    R = GRAYTR.Ray(Point3(1, 0, -0.5), Vector3(-1, 0, 0))
    dg, t, e = GRAYTR.shape_intersect(R, D)
    @test dg.p ≈ Point3(0.5, 0.0, -0.5)
    @test dg.n ≈ normalize(Normal3(1.0, 0.0, 1.0))
    @test t ≈ 0.5
    
    R = GRAYTR.Ray(Point3(0, 0, 1), Vector3(0, 0, -1), 10, Inf, 1)
    dg, t, e = GRAYTR.shape_intersect(R, D)
    @test dg == nothing
    @test isnan(t)
    @test isnan(e)
    
    R = GRAYTR.Ray(Point3(2.0, 0, 0.2), Vector3(-1, 0, 0))
    dg, t, e = GRAYTR.shape_intersect(R, D)
    @test dg == nothing
    @test isnan(t)
    @test isnan(e)
    
end # testset "intersections"



@testset "Create from point to point" begin

    cone = GRAYTR.cone_between_points(1, Point3(0,0,0), 1.0, Point3(0,0,1), 0.5)
    @test cone.h ≈ 0.5
    @test !any(isnan, cone.obj_to_world.M)
    @test !any(isnan, cone.world_to_obj.M)
    @test !any(isnan, cone.obj_to_world.MInv)
    @test !any(isnan, cone.world_to_obj.MInv)
    
    cone = GRAYTR.cone_between_points(1, Point3(0,0,0), 0.5, Point3(0,0,1), 1.0)
    @test cone.h ≈ 0.5
    @test !any(isnan, cone.obj_to_world.M)
    @test !any(isnan, cone.world_to_obj.M)
    @test !any(isnan, cone.obj_to_world.MInv)
    @test !any(isnan, cone.world_to_obj.MInv)
    
    cone = GRAYTR.cone_between_points(1, Point3(0,0,0), 1.0, Point3(0,0,-1), 0.5)
    @test cone.h ≈ 0.5
    @test !any(isnan, cone.obj_to_world.M)
    @test !any(isnan, cone.world_to_obj.M)
    @test !any(isnan, cone.obj_to_world.MInv)
    @test !any(isnan, cone.world_to_obj.MInv)
    
    cone = GRAYTR.cone_between_points(1, Point3(0,0,0), 0.5, Point3(0,0,-1), 1.0)
    @test cone.h ≈ 0.5
    @test !any(isnan, cone.obj_to_world.M)
    @test !any(isnan, cone.world_to_obj.M)
    @test !any(isnan, cone.obj_to_world.MInv)
    @test !any(isnan, cone.world_to_obj.MInv)
    
    

    P1 = Point3(1, 0, 0)
    P2 = Point3(3, 0, 0)
    r1 = 2.0
    r2 = 1.0

    cone = GRAYTR.cone_between_points(1, P1, r1, P2, r2)
    @test cone.h ≈ 0.5
    
    # Intersections from Y direction
    R = GRAYTR.Ray(Point3(0.0, -5.0, 0.0), Vector3(0, 1, 0))
    @test !GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, -5.0, 0.0), Vector3(0, 1, 0))
    @test GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, -5.0, 1.0), Vector3(0, 1, 0))
    @test GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, -5.0, -1.0), Vector3(0, 1, 0))
    @test GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(4.0, -5.0, 0.0), Vector3(0, 1, 0))
    @test !GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, -5.0, 2.0), Vector3(0, 1, 0))
    @test !GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, -5.0, -2.0), Vector3(0, 1, 0))
    @test !GRAYTR.intersectP(R, cone)
    
    
    # Intersections from Z direction
    R = GRAYTR.Ray(Point3(0.0, 0.0, 5.0), Vector3(0, 0, -1))
    @test !GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, 0.0, 5.0), Vector3(0, 0, -1))
    @test GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, -1.0, 5.0), Vector3(0, 0, -1))
    @test GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, -1.0, 5.0), Vector3(0, 0, -1))
    @test GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(4.0, 0.0, 5.0), Vector3(0, 0, -1))
    @test !GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, 2.0, 5.0), Vector3(0, 0, -1))
    @test !GRAYTR.intersectP(R, cone)
    
    R = GRAYTR.Ray(Point3(2.0, -2.0, 5.0), Vector3(0, 0, -1))
    @test !GRAYTR.intersectP(R, cone)
    
    
end




end # main Cone testset
