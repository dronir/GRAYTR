
@testset "Sphere" begin

S = Sphere()
@test S.id == 1
@test S.radius ≈ 1.0
@test GRAYTR.can_intersect(S)

S = Sphere(translation(-1, 0, 0))
S = Sphere(2.0, translation(-1, 0, 0))
S = Sphere(1, 2.0, translation(-1, 0, 0))
S = Sphere(1, 2.0, false, translation(-1, 0, 0))

@test GRAYTR.area(S) ≈ 16π

@test GRAYTR.obj_bounds(S) == GRAYTR.BoundingBox(Point3(-2,-2,-2), Point3(2,2,2))
@test GRAYTR.world_bounds(S) == GRAYTR.BoundingBox(Point3(-3,-2,-2), Point3(1,2,2))

S = Sphere()
T = scaling(2.0, 1.0, 1.0)
S = T(S)
@test GRAYTR.world_bounds(S) == GRAYTR.BoundingBox(Point3(-2,-1,-1), Point3(2,1,1))


end
