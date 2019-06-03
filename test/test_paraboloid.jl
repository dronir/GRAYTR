@testset "Paraboloid" begin

P = Paraboloid()

@test P.id == 1
@test P.h0 == 0.0
@test P.h1 == 1.0
@test P.radius == 1.0

# Create with transformation
P = Paraboloid(GRAYTR.scaling(2.0, 1.0, 1.0))
@test GRAYTR.obj_bounds(P) == GRAYTR.BoundingBox(Point3(-1, -1, 0), Point3(1, 1, 1))
@test GRAYTR.world_bounds(P) == GRAYTR.BoundingBox(Point3(-2, -1, 0), Point3(2, 1, 1))

# Transform after creation
T = GRAYTR.scaling(1.0, 2.0, 1.0)
P = T(Paraboloid())
@test GRAYTR.obj_bounds(P) == GRAYTR.BoundingBox(Point3(-1, -1, 0), Point3(1, 1, 1))
@test GRAYTR.world_bounds(P) == GRAYTR.BoundingBox(Point3(-1, -2, 0), Point3(1, 2, 1))


# intersectP
D = Paraboloid(translation(0, 0, -1))
R = GRAYTR.Ray(Point3(0, 0, 1), Vector3(0, 0, -1))
@test GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0, 0, 1), Vector3(0, 0, -1), 10, Inf, 1)
@test !GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0, 2, 1), Vector3(0, 0, -1))
@test !GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0.9*sqrt(2), 0.9*sqrt(2), 1), Vector3(0, 0, -1))
@test !GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0, 0, -2), Vector3(0, 0, 1))
@test GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0, 0, 0), Vector3(0, 1, 0))
@test GRAYTR.intersectP(R, D)


# intersect

R = GRAYTR.Ray(Point3(0, 0, 1), Vector3(0, 0, -1))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg.p ≈ Point3(0, 0, -1)
@test dg.n ≈ Normal3(0, 0, 1)
@test t ≈ 2.0

R = GRAYTR.Ray(Point3(0, 0, -3), Vector3(0, 0, 1))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg.p ≈ Point3(0, 0, -1)
@test dg.n ≈ Normal3(0, 0, -1)
@test t ≈ 2.0

R = GRAYTR.Ray(Point3(0.9*sqrt(2), 0.9*sqrt(2), 1), Vector3(0, 0, -1))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg == nothing
@test isnan(t)
@test isnan(e)

R = GRAYTR.Ray(Point3(1, 0, -0.5), Vector3(-1, 0, 0))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg.p ≈ Point3(1/sqrt(2), 0.0, -0.5)
@test dg.n ≈ Normal3(sqrt(2), 0.0, -1.0) / sqrt(3)
@test t ≈ 1 - 1/sqrt(2)

R = GRAYTR.Ray(Point3(0, 0, 1), Vector3(0, 0, -1), 10, Inf, 1)
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg == nothing
@test isnan(t)
@test isnan(e)


end
