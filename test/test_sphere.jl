
@testset "Sphere" begin

S = Sphere()
@test S.id == 1
@test S.radius ≈ 1.0
@test GRAYTR.can_intersect(S)

S = Sphere(translation(-1, 0, 0))
S = Sphere(2.0, π/2, translation(-1, 0, 0))
S = Sphere(1, 2.0, translation(-1, 0, 0))
S = Sphere(1, 2.0, π/2, translation(-1, 0, 0))

@test GRAYTR.area(S) ≈ 16π

@test GRAYTR.obj_bounds(S) == GRAYTR.BoundingBox(Point3(-2,-2,-2), Point3(2,2,2))
@test GRAYTR.world_bounds(S) == GRAYTR.BoundingBox(Point3(-3,-2,-2), Point3(1,2,2))

S = Sphere()
T = scaling(2.0, 1.0, 1.0)
S = T(S)
@test GRAYTR.world_bounds(S) == GRAYTR.BoundingBox(Point3(-2,-1,-1), Point3(2,1,1))


@testset "Full sphere intersections" begin



S = Sphere(1.0, π, Transformation())

R = Ray(Point3(2, 0, 0), Vector3(-1, 0, 0))
@test GRAYTR.intersectP(R, S)
R = Ray(Point3(-2, 0, 0), Vector3(1, 0, 0))
@test GRAYTR.intersectP(R, S)
R = Ray(Point3(2, 2, 0), Vector3(-1, 0, 0))
@test !GRAYTR.intersectP(R, S)
R = Ray(Point3(-2, 2, 0), Vector3(1, 0, 0))
@test !GRAYTR.intersectP(R, S)

R = Ray(Point3(0, 2, 0), Vector3(0, -1, 0))
@test GRAYTR.intersectP(R, S)
R = Ray(Point3(0, -2, 0), Vector3(0, 1, 0))
@test GRAYTR.intersectP(R, S)
R = Ray(Point3(0, 2, 2), Vector3(0, -1, 0))
@test !GRAYTR.intersectP(R, S)
R = Ray(Point3(0, -2, 2), Vector3(0, 1, 0))
@test !GRAYTR.intersectP(R, S)

R = Ray(Point3(0, 0, 2), Vector3(0, 0, -1))
@test GRAYTR.intersectP(R, S)
R = Ray(Point3(0, 0, -2), Vector3(0, 0, 1))
@test GRAYTR.intersectP(R, S)
R = Ray(Point3(2, 0, 2), Vector3(0, 0, -1))
@test !GRAYTR.intersectP(R, S)
R = Ray(Point3(2, 0, -2), Vector3(0, 0, 1))
@test !GRAYTR.intersectP(R, S)


# Shape intersect on a whole sphere

R = Ray(Point3(2, 0, 0), Vector3(-1, 0, 0))
dg, a, b = GRAYTR.shape_intersect(R, S)
@test dg.p ≈ Point3(1, 0, 0)
@test a ≈ 1

R = Ray(Point3(-2, 0, 0), Vector3(1, 0, 0))
dg, a, b = GRAYTR.shape_intersect(R, S)
@test dg.p ≈ Point3(-1, 0, 0)
@test a ≈ 1

R = Ray(Point3(2, 2, 0), Vector3(-1, 0, 0))
@test !GRAYTR.intersectP(R, S)
R = Ray(Point3(-2, 2, 0), Vector3(1, 0, 0))
dg, a, b = GRAYTR.shape_intersect(R, S)
@test dg == nothing
@test isnan(a)
@test isnan(b)

end

@testset "Hemisphere intersections" begin

# Shape intersection on a hemisphere

S = Sphere(1.0, π/2, Transformation())

R = Ray(Point3(0, 0, 2), Vector3(0, 0, -1))
dg, a, b = GRAYTR.shape_intersect(R, S)
@test dg.p ≈ Point3(0, 0, 1)
@test a ≈ 1

R = Ray(Point3(0, 0, -2), Vector3(0, 0, 1))
dg, a, b =  GRAYTR.shape_intersect(R, S)
@test dg.p ≈ Point3(0, 0, 1)
@test a ≈ 3

R = Ray(Point3(2, 0, -0.5), Vector3(-1, 0, 0))
dg, a, b =  GRAYTR.shape_intersect(R, S)
@test dg == nothing
@test isnan(a)
@test isnan(b)


end


end
