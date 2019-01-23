@testset "Cylinder" begin

C = GRAYTR.Cylinder(scaling(1, 1, 0.5))
@test C.id == 1
@test GRAYTR.can_intersect(C)
@test GRAYTR.area(C) ≈ 4π

@test GRAYTR.obj_bounds(C) == GRAYTR.BoundingBox(Point3(-1, -1, -1), Point3(1, 1, 1))
@test GRAYTR.world_bounds(C) == GRAYTR.BoundingBox(Point3(-1, -1, -0.5), Point3(1, 1, 0.5))

T = translation(0, 0, 0.5)
C2 = T(C)
@test GRAYTR.world_bounds(C2) == GRAYTR.BoundingBox(Point3(-1, -1, 0), Point3(1, 1, 1))

R = GRAYTR.Ray(Point3(5, 0, 0), Vector3(-1, 0, 0))
@test GRAYTR.intersectP(R, C)

R = GRAYTR.Ray(Point3(5, 4, 0), Vector3(-1, 0, 0))
@test !GRAYTR.intersectP(R, C)

R = GRAYTR.Ray(Point3(0, 0, 3), Vector3(0, 0, -1))
@test !GRAYTR.intersectP(R, C)

R = GRAYTR.Ray(Point3(5, 0, 0), Vector3(-1, 0, 0), 15.0, Inf, 1)
@test !GRAYTR.intersectP(R, C)


R = GRAYTR.Ray(Point3(5, 0, 0), Vector3(-1, 0, 0))
dg, t, e = GRAYTR.shape_intersect(R, C)
@test dg.p ≈ Point3(1, 0, 0)
@test t ≈ 4

R = GRAYTR.Ray(Point3(5, 4, 0), Vector3(-1, 0, 0))
dg, t, e = GRAYTR.shape_intersect(R, C)
@test dg == nothing
@test isnan(t)
@test isnan(e)

R = GRAYTR.Ray(Point3(0, 0, 3), Vector3(0, 0, -1))
dg, t, e = GRAYTR.shape_intersect(R, C)
@test dg == nothing
@test isnan(t)
@test isnan(e)


R = GRAYTR.Ray(Point3(5, 0, 0), Vector3(-1, 0, 0), 15.0, Inf, 1)
dg, t, e = GRAYTR.shape_intersect(R, C)
@test dg == nothing
@test isnan(t)
@test isnan(e)

end
