@testset "GeometricPrimitives" begin

S = GRAYTR.Sphere()
mat = GRAYTR.Lambert(NoLight())

G = GeometricPrimitive(S, mat, nothing, 1)
@test G.id == 1
@test G.shape == S
@test G.material == mat


shapes = [Sphere(), Disk()]
prims = apply_material(shapes, mat)
@test typeof(prims) == Array{GeometricPrimitive, 1}
@test length(prims) == 2
@test typeof(prims[1].shape) == Sphere
@test typeof(prims[2].shape) == Disk



end # main testset
