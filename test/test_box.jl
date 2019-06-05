@testset "Box" begin

box = Box()

@test box.id == 1
@test GRAYTR.can_intersect(box)

# Create with transformation
box = Box(GRAYTR.scaling(2.0, 1.0, 1.0))
@test GRAYTR.obj_bounds(box) == GRAYTR.BoundingBox(Point3(-1, -1, -1), Point3(1, 1, 1))
@test GRAYTR.world_bounds(box) == GRAYTR.BoundingBox(Point3(-2, -1, -1), Point3(2, 1, 1))

box = Box(2, GRAYTR.scaling(2.0, 1.0, 1.0))
@test GRAYTR.obj_bounds(box) == GRAYTR.BoundingBox(Point3(-1, -1, -1), Point3(1, 1, 1))
@test box.id == 2



# Transform after creation
T = GRAYTR.scaling(1.0, 2.0, 1.0)
box = T(Box())
@test GRAYTR.obj_bounds(box) == GRAYTR.BoundingBox(Point3(-1, -1, -1), Point3(1, 1, 1))
@test GRAYTR.world_bounds(box) == GRAYTR.BoundingBox(Point3(-1, -2, -1), Point3(1, 2, 1))


# intersectP
D = Box()

R = GRAYTR.Ray(Point3(0.5, 0, 2), Vector3(0, 0, -1))
@test GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0.5, 0, -2), Vector3(0, 0, 1))
@test GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0.5, 0, 2), Vector3(0, 0, -1), 10, Inf, 1)
@test !GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0.5, 0, -2), Vector3(0, 0, 1), 10, Inf, 1)
@test !GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0, 2, 2), Vector3(0, 0, -1))
@test !GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0, 2, -2), Vector3(0, 0, 1))
@test !GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0.5, 0, -2), Vector3(0, 0, 1))
@test GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0.5, 0, 2), Vector3(0, 0, -1))
@test GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0.0, 0, -0.2), Vector3(0, 1, 0))
@test GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0.0, 0, 0.0), Vector3(0, -1, 0))
@test GRAYTR.intersectP(R, D)


# intersect

R = GRAYTR.Ray(Point3(0.5, 0, 2), Vector3(0, 0, -1))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg.p ≈ Point3(0.5, 0, 1.0)
@test dg.n ≈ normalize(Normal3(0, 0, 1))
@test t ≈ 1.0

R = GRAYTR.Ray(Point3(0.5, 0, -3), Vector3(0, 0, 1))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg.p ≈ Point3(0.5, 0, -1.0)
@test dg.n ≈ normalize(Normal3(0.0, 0, -1.0))
@test t ≈ 2.0

R = GRAYTR.Ray(Point3(2, 0, -0.5), Vector3(-1, 0, 0))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg.p ≈ Point3(1.0, 0.0, -0.5)
@test dg.n ≈ normalize(Normal3(1.0, 0.0, 0.0))
@test t ≈ 1.0


# Ray along Y axis, hitting box

R = GRAYTR.Ray(Point3(0, 2, -0.5), Vector3(0, -1, 0))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg.p ≈ Point3(0.0, 1.0, -0.5)
@test dg.n ≈ normalize(Normal3(0.0, 1.0, 0.0))
@test t ≈ 1.0

R = GRAYTR.Ray(Point3(0, -2, -0.5), Vector3(0, 1, 0))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg.p ≈ Point3(0.0, -1.0, -0.5)
@test dg.n ≈ normalize(Normal3(0.0, -1.0, 0.0))
@test t ≈ 1.0

# Ray along Y axis, not hitting box

R = GRAYTR.Ray(Point3(0, 2, -1.5), Vector3(0, -1, 0))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg == nothing
@test isnan(t)
@test isnan(e)

R = GRAYTR.Ray(Point3(0, -2, -1.5), Vector3(0, 1, 0))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg == nothing
@test isnan(t)
@test isnan(e)


# Box not inside ray limits

R = GRAYTR.Ray(Point3(0, 0, 2), Vector3(0, 0, -1), 10, Inf, 1)
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg == nothing
@test isnan(t)
@test isnan(e)


end
