@testset "Triangle" begin

p1 = Point3(0)
p2 = Point3(1, 0, 0)
p3 = Point3(0, 1, 0)

Tri = GRAYTR.Triangle(p1, p2, p3)

@test GRAYTR.can_intersect(Tri)
@test GRAYTR.area(Tri) ≈ 0.5

@test GRAYTR.obj_bounds(Tri) == GRAYTR.BoundingBox(Point3(0), Point3(1, 1, 0))
@test GRAYTR.world_bounds(Tri) == GRAYTR.BoundingBox(Point3(0), Point3(1, 1, 0))

T = GRAYTR.translation(-0.5, -0.5, -1)
@test GRAYTR.world_bounds(T(Tri)) == GRAYTR.BoundingBox(Point3(-0.5, -0.5, -1), Point3(0.5, 0.5, -1))


R = GRAYTR.Ray(Point3(0.2, 0.2, 3), Vector3(0, 0, -1))
@test GRAYTR.intersectP(R, Tri)

R = GRAYTR.Ray(Point3(2, 1, 3), Vector3(0, 0, -1))
@test !GRAYTR.intersectP(R, Tri)

R = GRAYTR.Ray(Point3(2, 1, 3), Vector3(0, 0, 1))
@test !GRAYTR.intersectP(R, Tri)

R = GRAYTR.Ray(Point3(0.2, 0.2, 0), Vector3(1, 0, 0))
@test !GRAYTR.intersectP(R, Tri)

R = GRAYTR.Ray(Point3(0.2, 0.2, 0), Vector3(1, 0, 0), 5.0, Inf, 1)
@test !GRAYTR.intersectP(R, Tri)




R = GRAYTR.Ray(Point3(0.2, 0.2, 3), Vector3(0, 0, -1))
dg, t, e = GRAYTR.shape_intersect(R, Tri)
@test dg.n == Normal3(0, 0, 1)
@test t ≈ 3.0

R = GRAYTR.Ray(Point3(2, 1, 3), Vector3(0, 0, -1))
dg, t, e = GRAYTR.shape_intersect(R, Tri)
@test dg == nothing

R = GRAYTR.Ray(Point3(2, 1, 3), Vector3(0, 0, 1))
dg, t, e = GRAYTR.shape_intersect(R, Tri)
@test dg == nothing

R = GRAYTR.Ray(Point3(0.2, 0.2, 0), Vector3(1, 0, 0))
dg, t, e = GRAYTR.shape_intersect(R, Tri)
@test dg == nothing

R = GRAYTR.Ray(Point3(0.2, 0.2, 0), Vector3(1, 0, 0), 5.0, Inf, 1)
dg, t, e = GRAYTR.shape_intersect(R, Tri)
@test dg == nothing


end
