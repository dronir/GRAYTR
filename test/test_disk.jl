
@testset "Disk" begin

D = Disk()
@test D.id == 1
@test D.rmin ≈ 0
@test D.rmax ≈ 1
@test D.inverted == false
@test GRAYTR.area(D) ≈ π
@test GRAYTR.can_intersect(D)

# Create with transformation
D = Disk(GRAYTR.scaling(2.0, 1.0, 1.0))
@test GRAYTR.obj_bounds(D) == GRAYTR.BoundingBox(Point3(-1, -1, 0), Point3(1, 1, 0))
@test GRAYTR.world_bounds(D) == GRAYTR.BoundingBox(Point3(-2, -1, 0), Point3(2, 1, 0))

# Transform after creation
T = GRAYTR.scaling(1.0, 2.0, 1.0)
D = T(Disk())
@test GRAYTR.obj_bounds(D) == GRAYTR.BoundingBox(Point3(-1, -1, 0), Point3(1, 1, 0))
@test GRAYTR.world_bounds(D) == GRAYTR.BoundingBox(Point3(-1, -2, 0), Point3(1, 2, 0))

# intersectP
D = Disk(translation(0, 0, -1))
R = GRAYTR.Ray(Point3(0, 0, 1), Vector3(0, 0, -1))
@test GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0, 2, 1), Vector3(0, 0, -1))
@test !GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0.9*sqrt(2), 0.9*sqrt(2), 1), Vector3(0, 0, -1))
@test !GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0, 0, -1), Vector3(0, 0, 1))
@test GRAYTR.intersectP(R, D)

R = GRAYTR.Ray(Point3(0, 0, 0), Vector3(0, 1, 0))
@test !GRAYTR.intersectP(R, D)


# intersect

R = GRAYTR.Ray(Point3(0, 0, 1), Vector3(0, 0, -1))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg.p ≈ Point3(0, 0, -1)
@test dg.n ≈ Normal3(0, 0, 1)
@test t ≈ 2.0

R = GRAYTR.Ray(Point3(0.9*sqrt(2), 0.9*sqrt(2), 1), Vector3(0, 0, -1))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg == nothing
@test isnan(t)
@test isnan(e)

R = GRAYTR.Ray(Point3(1, 0, 0), Vector3(-1, 0, 0))
dg, t, e = GRAYTR.shape_intersect(R, D)
@test dg == nothing
@test isnan(t)
@test isnan(e)


end

